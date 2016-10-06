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

//________________________________________________________________
AliMultSelectionCuts::AliMultSelectionCuts() :
  TNamed(), fESD(0), fVzCut(10), fErrorCode(299),
    fEvSel_Trig_kMB (kTRUE),
    fEvSel_INELgtZERO (kTRUE),
    fEvSel_TrackletsVsClusters (kTRUE),
    fEvSel_RejectPileupInMultBins (kTRUE),
    fEvSel_CheckConsistencySPDandTrackVertices (kTRUE),
    fEvSel_NonZeroNContribs( kFALSE ),
    fEvSel_IsNotAsymmetricInVZERO( kFALSE ),
    fEvSel_IsNotIncompleteDAQ(kFALSE),
    fEvSel_HasGoodVertex2016(kFALSE)
{
  // Constructor
  
}
//________________________________________________________________
AliMultSelectionCuts::AliMultSelectionCuts(const char * name, const char * title):
  TNamed(name,title), fESD(0), fVzCut(10), fErrorCode(299),
    fEvSel_Trig_kMB (kTRUE),
    fEvSel_INELgtZERO (kTRUE),
    fEvSel_TrackletsVsClusters (kTRUE),
    fEvSel_RejectPileupInMultBins (kTRUE),
    fEvSel_CheckConsistencySPDandTrackVertices (kTRUE),
    fEvSel_NonZeroNContribs( kFALSE ),
    fEvSel_IsNotAsymmetricInVZERO( kFALSE ),
    fEvSel_IsNotIncompleteDAQ ( kFALSE ),
    fEvSel_HasGoodVertex2016(kFALSE)
{
  // Constructor
  
}
//________________________________________________________________
AliMultSelectionCuts& AliMultSelectionCuts::operator=(const AliMultSelectionCuts& o)
{
    if (&o == this) return *this;
    SetName(o.GetName());
    SetTitle(o.GetTitle());
    fESD = 0;
    fVzCut = o.fVzCut;
    fErrorCode = o.fErrorCode;
    fEvSel_Trig_kMB = o.fEvSel_Trig_kMB;
    fEvSel_INELgtZERO = o.fEvSel_INELgtZERO;
    fEvSel_TrackletsVsClusters = o.fEvSel_TrackletsVsClusters;
    fEvSel_RejectPileupInMultBins = o.fEvSel_RejectPileupInMultBins;
    fEvSel_CheckConsistencySPDandTrackVertices = o.fEvSel_CheckConsistencySPDandTrackVertices;
    fEvSel_NonZeroNContribs = o.fEvSel_NonZeroNContribs;
    fEvSel_IsNotAsymmetricInVZERO = o.fEvSel_IsNotAsymmetricInVZERO;
    fEvSel_IsNotIncompleteDAQ = o.fEvSel_IsNotIncompleteDAQ;
    fEvSel_HasGoodVertex2016 = o.fEvSel_HasGoodVertex2016;
    return *this;
}
//________________________________________________________________
AliMultSelectionCuts::~AliMultSelectionCuts(){
  // destructor
  
}

void AliMultSelectionCuts::Print(Option_t *option) const {
  Printf("Value of Cuts:");
  Printf(" Vertex Z position.....................: [%f]", fVzCut);
  Printf(" Physics Selection.....................: [%i]", fEvSel_Trig_kMB);
  Printf(" INEL > 0..............................: [%i]", fEvSel_INELgtZERO);
  Printf(" Tracklets vs Clusters.................: [%i]", fEvSel_TrackletsVsClusters);
  Printf(" Reject Pileup SPD (mult bins).........: [%i]", fEvSel_RejectPileupInMultBins);
  Printf(" SPD and Track vertex consistency......: [%i]", fEvSel_CheckConsistencySPDandTrackVertices);
  Printf(" Non Zero NContribs to PV..............: [%i]", fEvSel_NonZeroNContribs);
  Printf(" Reject Asymmetric in VZERO events.....: [%i]", fEvSel_IsNotAsymmetricInVZERO);
  Printf(" Reject IsIncompleteDAQ................: [%i]", fEvSel_IsNotIncompleteDAQ);
  Printf(" Reject events if vtx not good (2016)..: [%i]", fEvSel_HasGoodVertex2016);
}
//________________________________________________________________
/* Deprecated
 
 Bool_t AliMultSelectionCuts::IsEventSelected(AliESDEvent * esd) {
 //Possibly deprecated -> break down done inside AliMultSelectionTask
 fESD = esd;
 Bool_t returnValue = kTRUE;
 if(!IsMinBias())       returnValue = kFALSE;
 if(!IsVzCutSelected()) returnValue = kFALSE;
 if(!IsINELgtZERO())    returnValue = kFALSE;
 if(!HasNoInconsistentSPDandTrackVertices(esd)) returnValue = kFALSE;
 if( fESD -> IsPileupFromSPDInMultBins() ) returnValue = kFALSE;
 return returnValue;
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

Bool_t AliMultSelectionCuts::HasNoInconsistentSPDandTrackVertices(AliESDEvent * lESD)
{
    Bool_t returnValue = kTRUE;
    const AliESDVertex *lPrimaryVtxSPD    = NULL;
    const AliESDVertex *lPrimaryVtxTracks = NULL;
    
    lPrimaryVtxSPD    = lESD->GetPrimaryVertexSPD   ();
    lPrimaryVtxTracks = lESD->GetPrimaryVertexTracks();
    
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
} */
