//
// Class AliRsnCutEventUtils
//
// This cut implementation checks the quality of event primary vertex.
// Some functions currently work only with ESD events (not AOD).
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliRsnCutEventUtils.h"
#include "AliAnalysisUtils.h"
#include "AliESDtrackCuts.h"
#include "AliMultSelection.h"
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
  fCheckPileUpMultBins(kFALSE),
  fUseVertexSelection2013pA(kFALSE),
  fUseVertexSelection2013pAspectra(kFALSE),
  fMaxVtxZ(10.0),
  fFilterNSDeventsDPMJETpA2013(kFALSE),
  fCheckIncompleteDAQ(kFALSE),
  fCheckPastFuture(kFALSE),
  fPastFutureFirstBit(79),
  fPastFutureNBits(11),
  fCheckSPDClusterVsTrackletBG(kFALSE),
  fASPDCvsTCut(-9999.0),
  fBSPDCvsTCut(-9999.0),
  fCheckInelGt0SPDtracklets(kFALSE),
  fCheckAcceptedMultSelection(kFALSE),
  fUtils(0x0)
{
  //
  // Main constructor.
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
  fCheckPileUpMultBins(copy.fCheckPileUpMultBins),
  fUseVertexSelection2013pA(copy.fUseVertexSelection2013pA),
  fUseVertexSelection2013pAspectra(copy.fUseVertexSelection2013pAspectra),
  fMaxVtxZ(copy.fMaxVtxZ),
  fFilterNSDeventsDPMJETpA2013(copy.fFilterNSDeventsDPMJETpA2013),
  fCheckIncompleteDAQ(copy.fCheckIncompleteDAQ),
  fCheckPastFuture(copy.fCheckPastFuture),
  fPastFutureFirstBit(copy.fPastFutureFirstBit),
  fPastFutureNBits(copy.fPastFutureNBits),
  fCheckSPDClusterVsTrackletBG(copy.fCheckSPDClusterVsTrackletBG),
  fASPDCvsTCut(copy.fASPDCvsTCut),
  fBSPDCvsTCut(copy.fBSPDCvsTCut),
  fCheckInelGt0SPDtracklets(copy.fCheckInelGt0SPDtracklets),
  fCheckAcceptedMultSelection(copy.fCheckAcceptedMultSelection),
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
  fCheckPileUpMultBins=copy.fCheckPileUpMultBins;
  fUseVertexSelection2013pA=copy.fUseVertexSelection2013pA;
  fUseVertexSelection2013pAspectra=copy.fUseVertexSelection2013pAspectra;
  fMaxVtxZ=copy.fMaxVtxZ;
  fFilterNSDeventsDPMJETpA2013=copy.fFilterNSDeventsDPMJETpA2013;
  fCheckIncompleteDAQ=copy.fCheckIncompleteDAQ;
  fCheckPastFuture=copy.fCheckPastFuture;
  fPastFutureFirstBit=copy.fPastFutureFirstBit;
  fPastFutureNBits=copy.fPastFutureNBits;
  fCheckSPDClusterVsTrackletBG=copy.fCheckSPDClusterVsTrackletBG;
  fASPDCvsTCut=copy.fASPDCvsTCut;
  fBSPDCvsTCut=copy.fBSPDCvsTCut;
  fCheckInelGt0SPDtracklets=copy.fCheckInelGt0SPDtracklets;
  fCheckAcceptedMultSelection=copy.fCheckAcceptedMultSelection;
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
   if(fASPDCvsTCut>-9990.0) fUtils->SetASPDCvsTCut(fASPDCvsTCut);
   if(fBSPDCvsTCut>-9990.0) fUtils->SetBSPDCvsTCut(fBSPDCvsTCut);

   Bool_t accept = kTRUE;
   //remove first event in chunk 
   if ((fIsRmFirstEvInChunck) && (fUtils->IsFirstEventInChunk(vevt))) accept = kFALSE;

   //apply vertex selection - for 2013 pPb data
   if ((fUseVertexSelection2013pAspectra) && (!IsVertexSelected2013pAIDspectra(vevt))) accept = kFALSE;
   
   //apply vertex selection - for 2013 pPb data
   if((fUseVertexSelection2013pA) && (!fUtils->IsVertexSelected2013pA(vevt))) accept = kFALSE;
     
   // pile-up check
   if ((fCheckPileUppA2013) && (fUtils->IsPileUpEvent(vevt))) accept = kFALSE;
   if (fCheckPileUpMultBins && IsPileUpMultBins()) accept = kFALSE;

   //check if DAQ is incomplete
   if(fCheckIncompleteDAQ && IsIncompleteDAQ()) return kFALSE;

   //past/future protection
   if(fCheckPastFuture && FailsPastFuture()) return kFALSE;

   //check correlation between number of SPD clusters and number of tracklets
   if(fCheckSPDClusterVsTrackletBG && IsSPDClusterVsTrackletBG(vevt)) return kFALSE;

   //select INEL>0 based on SPD tracklets
   if(fCheckInelGt0SPDtracklets && !IsInelGt0SPDtracklets()) return kFALSE;

   if(fCheckAcceptedMultSelection && !IsAcceptedMultSelection()) return kFALSE;
   
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
       AliInfo("fFilterNSDeventsDPMJETpA2013 rejected a non-NSD event!");
       accept = kFALSE;
     }
     else  accept = kTRUE;   
   }
   
   return accept;
}


Bool_t AliRsnCutEventUtils::IsVertexSelected2013pAIDspectra(AliVEvent *event)
{
  // Check whether the event is of AOD or ESD type
  AliAODEvent *aod = 0x0;
  AliESDEvent *esd = 0x0;
  aod = dynamic_cast<AliAODEvent*>(event);
  esd = dynamic_cast<AliESDEvent*>(event);

  Bool_t fisAOD = kFALSE;
  if(aod) { 
    fisAOD = kTRUE; 
  } else {
    fisAOD = kFALSE;
    if(!esd) {
      AliFatal("Event is neither of AOD nor ESD type");
      return kFALSE;
    }
  }
  
  //Roberto's PV selection criteria, implemented 17th April 2013 on ESDs
  //fbellini modifications for usage with both ESDs and AODs
   /* vertex selection */
  
  Bool_t accept = kFALSE; 
  const AliVVertex *vertex = fisAOD ? 
    dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertex()) : 
    dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertex()) ;
  
  if (vertex->GetNContributors() < 1) {
    vertex = fisAOD ? 
      dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertexSPD()) : 
      dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertexSPD()) ;
    if (vertex->GetNContributors() < 1) accept = kFALSE;
    else accept = kTRUE;
    Double_t cov[6]={0};
    vertex->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vertex->IsFromVertexerZ() && (zRes>0.25)) accept = kFALSE;
  }
  else accept = kTRUE;    
  
  Float_t zvtx = vertex->GetZ();
  if (TMath::Abs(zvtx) > fMaxVtxZ) accept = kFALSE;
  return accept; 
}


Bool_t AliRsnCutEventUtils::IsPileUpMultBins(){
  if(!fEvent) return kFALSE;
  AliAODEvent* aodEvt=0;
  Bool_t isAOD=fEvent->IsAOD();
  if(isAOD) return kFALSE;//not yet available for AODs

  AliESDEvent* esdEvt=0;
  Bool_t isESD=fEvent->IsESD();
  if(isESD) esdEvt=dynamic_cast<AliESDEvent *>(fEvent->GetRef());

  else if(isESD && esdEvt && esdEvt->IsPileupFromSPDInMultBins()) return kTRUE;
  return kFALSE;
}


Bool_t AliRsnCutEventUtils::IsIncompleteDAQ(){
  if(!fEvent) return kFALSE;
  AliAODEvent* aodEvt=0;
  Bool_t isAOD=fEvent->IsAOD();
  if(isAOD) aodEvt=dynamic_cast<AliAODEvent *>(fEvent->GetRef());

  AliESDEvent* esdEvt=0;
  Bool_t isESD=fEvent->IsESD();
  if(isESD) esdEvt=dynamic_cast<AliESDEvent *>(fEvent->GetRef());

  if(isAOD && aodEvt && aodEvt->IsIncompleteDAQ()) return kTRUE;
  else if(isESD && esdEvt && esdEvt->IsIncompleteDAQ()) return kTRUE;
  return kFALSE;
}


Bool_t AliRsnCutEventUtils::FailsPastFuture(){
  if(!fEvent) return kFALSE;//this function has not been tested

  AliESDEvent* esdEvt=0;
  Bool_t isESD=fEvent->IsESD();
  if(isESD) esdEvt=(AliESDEvent*)(fEvent->GetRef());
  else return kFALSE;
  if(!esdEvt) return kFALSE;

  TBits b=esdEvt->GetHeader()->GetIRInt1InteractionMap();
  for(Int_t i=fPastFutureFirstBit;i<fPastFutureFirstBit+fPastFutureNBits;i++) if(b.TestBitNumber(i)) return kTRUE;
  return kFALSE;
}


Bool_t AliRsnCutEventUtils::IsSPDClusterVsTrackletBG(AliVEvent* vevt){
  if(!fEvent || !fUtils) return kFALSE;
  if(!vevt) vevt = dynamic_cast<AliVEvent *>(fEvent->GetRef());
  if(!vevt) return kFALSE;
  return fUtils->IsSPDClusterVsTrackletBG(vevt);
}


Bool_t AliRsnCutEventUtils::IsInelGt0SPDtracklets(){
  if(!fEvent) return kFALSE;
  AliAODEvent* aodEvt=0;
  Bool_t isAOD=fEvent->IsAOD();
  if(isAOD) aodEvt=dynamic_cast<AliAODEvent *>(fEvent->GetRef());

  AliESDEvent* esdEvt=0;
  Bool_t isESD=fEvent->IsESD();
  if(isESD) esdEvt=dynamic_cast<AliESDEvent *>(fEvent->GetRef());

  if(isAOD){
    if(aodEvt && IsInelGt0SPDtracklets(aodEvt)) return kTRUE;
    else return kFALSE;
  }else if(isESD){
    if(esdEvt && IsInelGt0SPDtracklets(esdEvt)) return kTRUE;
    else return kFALSE;
  }
  return kFALSE;
}


Bool_t AliRsnCutEventUtils::IsInelGt0SPDtracklets(AliAODEvent* aodEvt){
  return kTRUE;//not yet implemented for AODs
}


Bool_t AliRsnCutEventUtils::IsInelGt0SPDtracklets(AliESDEvent* esdEvt){
  if(!esdEvt) return kTRUE;
  Int_t n=AliESDtrackCuts::GetReferenceMultiplicity(esdEvt,AliESDtrackCuts::kTracklets,1.,0.);
  if(n<1) return kFALSE;
  return kTRUE;
}


Bool_t AliRsnCutEventUtils::IsAcceptedMultSelection(){
  AliMultSelection *MultSelection=0;

  if(!fEvent) return kFALSE;
  AliAODEvent* aodEvt=0;
  Bool_t isAOD=fEvent->IsAOD();
  if(isAOD) aodEvt=dynamic_cast<AliAODEvent *>(fEvent->GetRef());
  if(isAOD && aodEvt){
    MultSelection=(AliMultSelection*) aodEvt->FindListObject("MultSelection");
    if(!MultSelection) return kTRUE;
    if(MultSelection->IsEventSelected()) return kTRUE;
    return kFALSE;
  }

  AliESDEvent* esdEvt=0;
  Bool_t isESD=fEvent->IsESD();
  if(isESD) esdEvt=dynamic_cast<AliESDEvent *>(fEvent->GetRef());
  if(isESD && esdEvt){
    MultSelection=(AliMultSelection*) esdEvt->FindListObject("MultSelection");
    if(!MultSelection) return kTRUE;
    if(MultSelection->IsEventSelected()) return kTRUE;
    return kFALSE;
  }

  return kFALSE;
}
