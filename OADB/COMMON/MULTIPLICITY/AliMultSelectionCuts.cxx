/**********************************************
 *
 *   Class meant to manage various 
 *   event selection switches, configurations
 *   and procedures
 *
 **********************************************/

#include "AliMultSelectionCuts.h"

#if !defined (__CINT__) || (defined(__MAKECINT__))
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#endif

ClassImp(AliMultSelectionCuts);

AliMultSelectionCuts::AliMultSelectionCuts() :
  TNamed(), fESD(0), fVzCut(10), fErrorCode(0),
    fEvSel_Trig_kMB (kTRUE),
    fEvSel_INELgtZERO (kTRUE),
    fEvSel_TrackletsVsClusters (kTRUE),
    fEvSel_RejectPileupInMultBins (kTRUE),
    fEvSel_CheckConsistencySPDandTrackVertices (kTRUE)
{
  // Constructor
  
}
AliMultSelectionCuts::AliMultSelectionCuts(const char * name, const char * title):
  TNamed(name,title), fESD(0), fVzCut(10), fErrorCode(0),
    fEvSel_Trig_kMB (kTRUE),
    fEvSel_INELgtZERO (kTRUE),
    fEvSel_TrackletsVsClusters (kTRUE),
    fEvSel_RejectPileupInMultBins (kTRUE),
    fEvSel_CheckConsistencySPDandTrackVertices (kTRUE) 
{
  // Constructor
  
}
AliMultSelectionCuts::~AliMultSelectionCuts(){
  // destructor
  
}
Bool_t AliMultSelectionCuts::IsEventSelected(AliESDEvent * esd) {
  fESD = esd;
  Bool_t returnValue = kTRUE;
    if(!IsMinBias())       returnValue = kFALSE;
    if(!IsVzCutSelected()) returnValue = kFALSE;
    if(!IsINELgtZERO())    returnValue = kFALSE;
    if(!HasNoInconsistentSPDandTrackVertices()) returnValue = kFALSE;
    return returnValue;
}

void AliMultSelectionCuts::Print(Option_t *option) const {
  Printf("Value of Cuts:");
  Printf(" - Vz: [%f]", fVzCut);
}

Bool_t AliMultSelectionCuts::IsMinBias(){
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
    return isSelected;
}


Bool_t AliMultSelectionCuts::IsVzCutSelected(){
    const AliESDVertex *lPrimaryBestESDVtx     = fESD->GetPrimaryVertex();
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    Bool_t returnValue = kTRUE;
    if (TMath::Abs( lBestPrimaryVtxPos[2] ) > fVzCut ) returnValue = kFALSE;
    return returnValue;
}

Bool_t AliMultSelectionCuts::IsINELgtZERO(){
    Bool_t returnValue = kFALSE;
    if ( AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0) >= 1 ) returnValue = kTRUE;
    return returnValue;
}

Bool_t AliMultSelectionCuts::HasNoInconsistentSPDandTrackVertices()
{
    Bool_t returnValue = kTRUE;
    const AliESDVertex *lPrimaryVtxSPD    = NULL;
    const AliESDVertex *lPrimaryVtxTracks = NULL;
    
    lPrimaryVtxSPD    = fESD->GetPrimaryVertexSPD   ();
    lPrimaryVtxTracks = fESD->GetPrimaryVertexTracks();
    
    //Only continue if track vertex defined
    if( lPrimaryVtxTracks->GetStatus() && lPrimaryVtxSPD->GetStatus() ){
        //Copy-paste from refmult estimator
        // TODO value of displacement to be studied
        const Float_t maxDisplacement = 0.5;
        //check for displaced vertices
        Double_t displacement = TMath::Abs(lPrimaryVtxSPD->GetZ() - lPrimaryVtxTracks->GetZ());
        if (displacement > maxDisplacement) returnValue = kFALSE;
    }
    return returnValue;
}
