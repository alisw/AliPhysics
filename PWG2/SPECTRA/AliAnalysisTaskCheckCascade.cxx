
/**************************************************************************
 *  Authors : Antonin Maire, Boris Hippolyte                              *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//                 AliAnalysisTaskCheckCascade class
//            (AliAnalysisTaskCheckCascade)
//            This task has four roles :
//              1. QAing the Cascades from ESD and AOD
//                 Origin:  AliAnalysisTaskESDCheckV0 by B.H. Nov2007, hippolyt@in2p3.fr
//              2. Prepare the plots which stand as raw material for yield extraction (wi/wo PID)
//              3. Supply an AliCFContainer meant to define the optimised topological selections
//              4. Rough azimuthal correlation study (Eta, Phi)
//            Adapted to Cascade : A.Maire Mar2008, antonin.maire@ires.in2p3.fr
//            Modified :           A.Maire Nov2010, antonin.maire@ires.in2p3.fr
//-----------------------------------------------------------------



class TTree;
class TParticle;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"


#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
//     #include "AliV0vertexer.h"
//     #include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"

#include "AliAnalysisTaskCheckCascade.h"

ClassImp(AliAnalysisTaskCheckCascade)



//________________________________________________________________________
AliAnalysisTaskCheckCascade::AliAnalysisTaskCheckCascade() 
  : AliAnalysisTaskSE(), fAnalysisType("ESD"), fTriggerMaskType("kMB"), fCollidingSystems(0), fESDpid(0), fESDtrackCuts(0), /*fPaveTextBookKeeping(0),*/
    fkRerunV0CascVertexers         (0),
    fkQualityCutZprimVtxPos        (kTRUE),
    fkRejectEventPileUp            (kTRUE),
    fkQualityCutNoTPConlyPrimVtx   (kTRUE),
    fkQualityCutTPCrefit           (kTRUE),
    fkQualityCut80TPCcls           (kTRUE),
    fkIsDataRecoWith1PadTPCCluster (kTRUE),
    fkExtraSelections              (0),
    fAngularCorrelationType        ("TrigLeadingTrck-AssoCasc"),

    	// - Cascade part initialisation
    fListHistCascade(0),
    fHistCascadeMultiplicityBeforeTrigSel(0),
    fHistCascadeMultiplicityForTrigEvt(0), fHistTrackMultiplicityForTrigEvt(0), fHistTPCrefitTrackMultiplicityForTrigEvt(0), fHistPrimaryTrackMultiplicityForTrigEvt(0),
    fHistCascadeMultiplicityForTrigEvtAndZprimVtx(0), fHistCascadeMultiplicityForTrigEvtNonPiledUpAndZprimVtx(0),
    fHistCascadeMultiplicityForSelEvt(0),
    fHistPosBestPrimaryVtxXForSelEvt(0), fHistPosBestPrimaryVtxYForSelEvt(0), fHistPosBestPrimaryVtxZForSelEvt(0),
    fHistTPCrefitTrackMultiplicityForCascadeEvt(0), fHistPrimaryTrackMultiplicityForCascadeEvt(0),
    fHistPosV0TPCClusters(0), fHistNegV0TPCClusters(0), fHistBachTPCClusters(0),
    fHistVtxStatus(0),

    fHistPosTrkgPrimaryVtxXForCascadeEvt(0), fHistPosTrkgPrimaryVtxYForCascadeEvt(0), fHistPosTrkgPrimaryVtxZForCascadeEvt(0), fHistTrkgPrimaryVtxRadius(0),
    fHistPosBestPrimaryVtxXForCascadeEvt(0), fHistPosBestPrimaryVtxYForCascadeEvt(0), fHistPosBestPrimaryVtxZForCascadeEvt(0), fHistBestPrimaryVtxRadius(0),
    f2dHistTrkgPrimVtxVsBestPrimVtx(0),

    fHistEffMassXi(0),  fHistChi2Xi(0),  
    fHistDcaXiDaughters(0), fHistDcaBachToPrimVertex(0), fHistXiCosineOfPointingAngle(0), fHistXiRadius(0),

    fHistMassLambdaAsCascDghter(0),
    fHistV0Chi2Xi(0),
    fHistDcaV0DaughtersXi(0),
    fHistDcaV0ToPrimVertexXi(0), 
    fHistV0CosineOfPointingAngleXi(0),
    fHistV0RadiusXi(0),
    fHistDcaPosToPrimVertexXi(0), fHistDcaNegToPrimVertexXi(0), 

    fHistMassXiMinus(0), fHistMassXiPlus(0),
    fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
    fHistMassWithCombPIDXiMinus(0), fHistMassWithCombPIDXiPlus(0),
    fHistMassWithCombPIDOmegaMinus(0), fHistMassWithCombPIDOmegaPlus(0),

    fHistXiTransvMom(0),    fHistXiTotMom(0),
    fHistBachTransvMomXi(0),   fHistBachTotMomXi(0),

    fHistChargeXi(0),
    fHistV0toXiCosineOfPointingAngle(0),

    fHistRapXi(0), fHistRapOmega(0), fHistEtaXi(0),
    fHistThetaXi(0), fHistPhiXi(0),

    fHistcTauXiMinus(0), fHistcTauXiPlus(0), fHistcTauOmegaMinus(0), fHistcTauOmegaPlus(0),

    f2dHistArmenteros(0),			
    f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
    f2dHistEffMassLambdaVsEffMassXiPlus(0), f2dHistEffMassXiVsEffMassOmegaPlus(0),
    f2dHistXiRadiusVsEffMassXiMinus(0), f2dHistXiRadiusVsEffMassXiPlus(0),
    f2dHistXiRadiusVsEffMassOmegaMinus(0), f2dHistXiRadiusVsEffMassOmegaPlus(0),
    
    f2dHistTPCdEdxOfCascDghters(0),
    
    f3dHistXiPtVsEffMassVsYXiMinus(0), f3dHistXiPtVsEffMassVsYXiPlus(0),
    f3dHistXiPtVsEffMassVsYOmegaMinus(0), f3dHistXiPtVsEffMassVsYOmegaPlus(0),
    
    fCFContCascadePIDXiMinus(0),
    fCFContCascadePIDXiPlus(0),
    fCFContCascadePIDOmegaMinus(0),
    fCFContCascadePIDOmegaPlus(0),
    fCFContCascadeCuts(0),
    
    fHnSpAngularCorrXiMinus(0), fHnSpAngularCorrXiPlus(0), 
    fHnSpAngularCorrOmegaMinus(0), fHnSpAngularCorrOmegaPlus(0)

{
  // Dummy Constructor
        for(Int_t iAlephIdx   = 0; iAlephIdx   < 5; iAlephIdx++   ) { fAlephParameters [iAlephIdx]    = -1.; }
        for(Int_t iV0selIdx   = 0; iV0selIdx   < 7; iV0selIdx++   ) { fV0Sels          [iV0selIdx   ] = -1.; }
        for(Int_t iCascSelIdx = 0; iCascSelIdx < 8; iCascSelIdx++ ) { fCascSels        [iCascSelIdx ] = -1.; }
}








//________________________________________________________________________
AliAnalysisTaskCheckCascade::AliAnalysisTaskCheckCascade(const char *name) 
  : AliAnalysisTaskSE(name), fAnalysisType("ESD"), fTriggerMaskType("kMB"), fCollidingSystems(0), fESDpid(0), fESDtrackCuts(0), /*fPaveTextBookKeeping(0),*/
    fkRerunV0CascVertexers         (0),
    fkQualityCutZprimVtxPos        (kTRUE),
    fkRejectEventPileUp            (kTRUE),
    fkQualityCutNoTPConlyPrimVtx   (kTRUE),
    fkQualityCutTPCrefit           (kTRUE),
    fkQualityCut80TPCcls           (kTRUE),
    fkIsDataRecoWith1PadTPCCluster (kTRUE),
    fkExtraSelections              (0),
    fAngularCorrelationType        ("TrigLeadingTrck-AssoCasc"),
     
    	// - Cascade part initialisation
    fListHistCascade(0),
    fHistCascadeMultiplicityBeforeTrigSel(0),
    fHistCascadeMultiplicityForTrigEvt(0), fHistTrackMultiplicityForTrigEvt(0), fHistTPCrefitTrackMultiplicityForTrigEvt(0), fHistPrimaryTrackMultiplicityForTrigEvt(0),
    fHistCascadeMultiplicityForTrigEvtAndZprimVtx(0), fHistCascadeMultiplicityForTrigEvtNonPiledUpAndZprimVtx(0),
    fHistCascadeMultiplicityForSelEvt(0),
    fHistPosBestPrimaryVtxXForSelEvt(0), fHistPosBestPrimaryVtxYForSelEvt(0), fHistPosBestPrimaryVtxZForSelEvt(0),
    fHistTPCrefitTrackMultiplicityForCascadeEvt(0), fHistPrimaryTrackMultiplicityForCascadeEvt(0),
    fHistPosV0TPCClusters(0), fHistNegV0TPCClusters(0), fHistBachTPCClusters(0),
    fHistVtxStatus(0),

    fHistPosTrkgPrimaryVtxXForCascadeEvt(0), fHistPosTrkgPrimaryVtxYForCascadeEvt(0), fHistPosTrkgPrimaryVtxZForCascadeEvt(0), fHistTrkgPrimaryVtxRadius(0),
    fHistPosBestPrimaryVtxXForCascadeEvt(0), fHistPosBestPrimaryVtxYForCascadeEvt(0), fHistPosBestPrimaryVtxZForCascadeEvt(0), fHistBestPrimaryVtxRadius(0),
    f2dHistTrkgPrimVtxVsBestPrimVtx(0),

    fHistEffMassXi(0),  fHistChi2Xi(0),  
    fHistDcaXiDaughters(0), fHistDcaBachToPrimVertex(0), fHistXiCosineOfPointingAngle(0), fHistXiRadius(0),

    fHistMassLambdaAsCascDghter(0),
    fHistV0Chi2Xi(0),
    fHistDcaV0DaughtersXi(0),
    fHistDcaV0ToPrimVertexXi(0), 
    fHistV0CosineOfPointingAngleXi(0),
    fHistV0RadiusXi(0),
    fHistDcaPosToPrimVertexXi(0), fHistDcaNegToPrimVertexXi(0), 

    fHistMassXiMinus(0), fHistMassXiPlus(0),
    fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
    fHistMassWithCombPIDXiMinus(0), fHistMassWithCombPIDXiPlus(0),
    fHistMassWithCombPIDOmegaMinus(0), fHistMassWithCombPIDOmegaPlus(0),

    fHistXiTransvMom(0),    fHistXiTotMom(0),
    fHistBachTransvMomXi(0),   fHistBachTotMomXi(0),

    fHistChargeXi(0),
    fHistV0toXiCosineOfPointingAngle(0),

    fHistRapXi(0), fHistRapOmega(0), fHistEtaXi(0),
    fHistThetaXi(0), fHistPhiXi(0),

    fHistcTauXiMinus(0), fHistcTauXiPlus(0), fHistcTauOmegaMinus(0), fHistcTauOmegaPlus(0),

    f2dHistArmenteros(0),			
    f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
    f2dHistEffMassLambdaVsEffMassXiPlus(0), f2dHistEffMassXiVsEffMassOmegaPlus(0),
    f2dHistXiRadiusVsEffMassXiMinus(0), f2dHistXiRadiusVsEffMassXiPlus(0),
    f2dHistXiRadiusVsEffMassOmegaMinus(0), f2dHistXiRadiusVsEffMassOmegaPlus(0),
    
    f2dHistTPCdEdxOfCascDghters(0),
    
    f3dHistXiPtVsEffMassVsYXiMinus(0), f3dHistXiPtVsEffMassVsYXiPlus(0),
    f3dHistXiPtVsEffMassVsYOmegaMinus(0), f3dHistXiPtVsEffMassVsYOmegaPlus(0),
    
    fCFContCascadePIDXiMinus(0),
    fCFContCascadePIDXiPlus(0),
    fCFContCascadePIDOmegaMinus(0),
    fCFContCascadePIDOmegaPlus(0),
    fCFContCascadeCuts(0),
    
    fHnSpAngularCorrXiMinus(0), fHnSpAngularCorrXiPlus(0), 
    fHnSpAngularCorrOmegaMinus(0), fHnSpAngularCorrOmegaPlus(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #1 writes into a TList container (cascade)
        
        for(Int_t iAlephIdx   = 0; iAlephIdx   < 5; iAlephIdx++   ) { fAlephParameters [iAlephIdx]    = -1.; }
        
        // Loose
        /*
//         fV0Sels[0] =  33.  ;  // max allowed chi2
//         fV0Sels[1] =   0.01;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
//         fV0Sels[2] =   0.01;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
//         fV0Sels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
//         fV0Sels[4] =   0.0 ;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
//         fV0Sels[5] =   0.2 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
//         fV0Sels[6] = 100.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
//         
//         fCascSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
//         fCascSels[1] =   0.02 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
//         fCascSels[2] =   0.008;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
//         fCascSels[3] =   0.01 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
//         fCascSels[4] =   0.5  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
//         fCascSels[5] =   0.98 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
//         fCascSels[6] =   0.2  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
//         fCascSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
        */
        
        
        // Hyper Loose "à la 900 GeV 2009 data", with lower cosine of pointing angle for Xi (0.95 down to 0.82) = 900 GeV paper
        
        fV0Sels[0] =  33.  ;  // max allowed chi2
        fV0Sels[1] =   0.001; // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
        fV0Sels[2] =   0.001; // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
        fV0Sels[3] =   5.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
        fV0Sels[4] =   0.0 ;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
        fV0Sels[5] =   0.1 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
        fV0Sels[6] = 100.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
        
        fCascSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
        fCascSels[1] =   0.001;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
        fCascSels[2] =   0.008;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
        fCascSels[3] =   0.001;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
        fCascSels[4] =   5.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
        fCascSels[5] =   0.82 ;  //FIXME min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
        fCascSels[6] =   0.1  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
        fCascSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
        
        
        //New default vtxR (http://alisoft.cern.ch/viewvc?view=rev&root=AliRoot&revision=40955, 5 May 2010)
        /*
        fV0Sels[0] =  33.  ;  // max allowed chi2
        fV0Sels[1] =   0.05;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
        fV0Sels[2] =   0.05;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
        fV0Sels[3] =   1.5 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
        fV0Sels[4] =   0.9 ;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
        fV0Sels[5] =   0.2 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
        fV0Sels[6] = 100.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
        
        fCascSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
        fCascSels[1] =   0.01 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
        fCascSels[2] =   0.008;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
        fCascSels[3] =   0.01 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
        fCascSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
        fCascSels[5] =   0.98 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
        fCascSels[6] =   0.2  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
        fCascSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
        */
        
        // Tight for Xi in p-p (http://alisoft.cern.ch/viewvc?view=rev&root=AliRoot&revision=40955, 5 May 2010)
        /*
        fV0Sels[0] =  33.  ;  // max allowed chi2
        fV0Sels[1] =   0.05;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
        fV0Sels[2] =   0.05;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
        fV0Sels[3] =   0.5 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
        fV0Sels[4] =   0.99 ; // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
        fV0Sels[5] =   0.2 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
        fV0Sels[6] = 100.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
        
        fCascSels[0] =  33.   ;   // max allowed chi2 (same as PDC07)
        fCascSels[1] =   0.05 ;   // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
        fCascSels[2] =   0.01;    // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
        fCascSels[3] =   0.035 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
        fCascSels[4] =   0.10  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
        fCascSels[5] =   0.9985 ; // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
        fCascSels[6] =   0.2  ;   // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
        fCascSels[7] = 100.   ;   // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
        */
        
        //NOTE
        // For PbPb coming data, take a look at Iouri's proposal
        // https://savannah.cern.ch/bugs/index.php?69877

        
  // Output slot #0 writes into a TList container (Cascade)
  DefineOutput(1, TList::Class());
  /*DefineOutput(2, TPaveText::Class());*/
}


AliAnalysisTaskCheckCascade::~AliAnalysisTaskCheckCascade()
{
  //
  // Destructor
  //

  // For all TH1, 2, 3 HnSparse and CFContainer are in the fListCascade TList.
  // They will be deleted when fListCascade is deleted by the TSelector dtor
  // Because of TList::SetOwner() ...
        
  if (fListHistCascade)         { delete fListHistCascade;     fListHistCascade = 0x0;    }
  if (fESDpid)                  { delete fESDpid;              fESDpid = 0x0;} // fESDpid is not stored into the TList
  if (fESDtrackCuts)            { delete fESDtrackCuts;        fESDtrackCuts = 0x0; }
  //if (fPaveTextBookKeeping)     { delete fPaveTextBookKeeping; fPaveTextBookKeeping = 0x0;} // fPaveTextBookKeeping is not strored in the TList
}



//________________________________________________________________________
void AliAnalysisTaskCheckCascade::UserCreateOutputObjects()
{
  // Create histograms
  // Called once



 fListHistCascade = new TList();
 fListHistCascade->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
 
 
if(! fESDpid){

  AliMCEventHandler *lmcEvtHandler  = dynamic_cast<AliMCEventHandler*>( (AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler() );
  
  if( !lmcEvtHandler ){ // !0x0 = real data or !1 = there is an MC handler available (useMC = kTRUE in AnalysisTrainNew), so = data from MC
        
        if(fkIsDataRecoWith1PadTPCCluster){
                // Parameters extracted for LHC10d Pass2 - A. Kalweit
                fAlephParameters[0] = 1.28949/50.; // = 0,0257898
                fAlephParameters[1] = 2.74095e+01;
                fAlephParameters[2] = TMath::Exp(-3.21763e+01);
                fAlephParameters[3] = 2.44026;
                fAlephParameters[4] = 6.58800;
        }
        else {
                // Reasonable parameters extracted for real p-p event (Dec 2009 - GSI Pass5) - A.Kalweit
                fAlephParameters[0] = 0.0283086;        // No extra-division to apply in SetBlochParam
                fAlephParameters[1] = 2.63394e+01;
                fAlephParameters[2] = 5.04114e-11;
                fAlephParameters[3] = 2.12543e+00;
                fAlephParameters[4] = 4.88663e+00; 
        }
        
        Printf("CheckCascade -      Check Aleph Param in case of REAL Data (fAlephParameters[3] = %f)  (To be compared with : 2.44026 for 1-pad-cluster prod. / 2.12543 otherwise)\n",  fAlephParameters[3]);

  }
  else { // MC reco
        if(fkIsDataRecoWith1PadTPCCluster){
                // Home made parameterization for LHC10f6a production = p+p 7 TeV
                fAlephParameters[0] = 0.04;
                fAlephParameters[1] = 17.5;
                fAlephParameters[2] = 3.4e-09;
                fAlephParameters[3] = 2.15;
                fAlephParameters[4] = 3.91720e+00; 
                
                // Home made parameterization for LHC10e13 production = p+p 900 GeV/c
                
        }
        else {
                // Reasonable parameters extracted for p-p simulation (LHC09a4) - A.Kalweit
                // fAlephParameters[0] = 4.23232575531564326e+00;//50*0.76176e-1; // do not forget to divide this value by 50 in SetBlochParam !
                // fAlephParameters[1] = 8.68482806165147636e+00;//10.632; 
                // fAlephParameters[2] = 1.34000000000000005e-05;//0.13279e-4;
                // fAlephParameters[3] = 2.30445734159456084e+00;//1.8631;
                // fAlephParameters[4] = 2.25624744086878559e+00;//1.9479;  
                
                // Reasonable parameters extracted for MC LHC09d10 event (Jan 2010) - A.Kalweit
                fAlephParameters[0] = 2.15898e+00/50.;
                fAlephParameters[1] = 1.75295e+01;
                fAlephParameters[2] = 3.40030e-09;
                fAlephParameters[3] = 1.96178e+00;
                fAlephParameters[4] = 3.91720e+00;
        }
        Printf("CheckCascade -      Check Aleph Param in case of MC Data (fAlephParameters[3] = %f) (To be compared with : 2.15 for 1-pad-cluster prod. / 1.96178 otherwise)\n",  fAlephParameters[3]);
  }

  fESDpid = new AliESDpid();
  fESDpid->GetTPCResponse().SetBetheBlochParameters( fAlephParameters[0],
                                                     fAlephParameters[1],
                                                     fAlephParameters[2],
                                                     fAlephParameters[3],
                                                     fAlephParameters[4] );
}
 
if(! fESDtrackCuts ){
      fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE); // Std definition of primary (see kTRUE argument) tracks for 2010
      fESDtrackCuts->SetEtaRange(-0.8,+0.8);
      fESDtrackCuts->SetPtRange(0.15, 1e10);
      Printf("CheckCascade - ESDtrackCuts set up to 2010 std ITS-TPC cuts...");
}
 
/* 
if( !fPaveTextBookKeeping){
        // FIXME : prepare a field with the AliRoot+root distrib ...

        fPaveTextBookKeeping = new TPaveText(0.1, 0.1, 0.9, 0.9,"NDC");
        fPaveTextBookKeeping->SetName("fPaveTextBookKeeping");
        fPaveTextBookKeeping->SetBorderSize(0);
        fPaveTextBookKeeping->SetTextAlign(12);
        fPaveTextBookKeeping->SetFillColor(kWhite);
        fPaveTextBookKeeping->SetTextFont(42);        // regular Arial or Helvetica,
        fPaveTextBookKeeping->SetTextColor(kBlue+4);
        
        fPaveTextBookKeeping->AddText( "Task CHECK CASCADE analysis" );
        fPaveTextBookKeeping->AddText("- - - - - - - - - - - ");
        fPaveTextBookKeeping->AddText( Form("AnalysisType : %s ", fAnalysisType.Data() ));
        if(!fCollidingSystems)  fPaveTextBookKeeping->AddText("Colliding system : p-p collisions ");
        else                    fPaveTextBookKeeping->AddText("Colliding system : A-A collisions ");

        fPaveTextBookKeeping->AddText("- - - - - - - - - - - ");
    
        if(fkRerunV0CascVertexers){
                fPaveTextBookKeeping->AddText("A.1. With V0 vertexer : ");
                fPaveTextBookKeeping->AddText( Form("  - V0 #chi^{2} _________________ <  %.3f ",               fV0Sels[0] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA(prim. Vtx/ 1^{st} daughter) ___ >  %.3f     cm ",  fV0Sels[1] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA(prim. Vtx/ 2^{nd} daughter) __  >  %.3f     cm",   fV0Sels[2] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA between V0 daughters ___ <  %.3f      cm",         fV0Sels[3] ));
                fPaveTextBookKeeping->AddText( Form("  - cos(V0 pointing angle) _______ >  %.3f ",              fV0Sels[4] ));
                fPaveTextBookKeeping->AddText( Form("  - R_{transv}(V0 decay) ________ >  %.3f             cm", fV0Sels[5] ));
                fPaveTextBookKeeping->AddText( Form("  - R_{transv}(V0 decay) ________ <  %.3f         cm",     fV0Sels[6] ));
                
                fPaveTextBookKeeping->AddText(" "); 
                
                fPaveTextBookKeeping->AddText("A.2. With Casc. vertexer : ");
                fPaveTextBookKeeping->AddText( Form("  - Casc. #chi^{2} ______________  <  %.3f ",                               fCascSels[0] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA(prim. Vtx/ V0) _________ >  %.3f    cm",                            fCascSels[1] ));
                fPaveTextBookKeeping->AddText( Form("  - | M_{#Lambda}(reco) - M_{#Lambda}(pdg) | _______ <  %.3f    GeV/c^{2}", fCascSels[2] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA(prim. Vtx/ Bach) _______ >  %.3f    cm",                            fCascSels[3] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA between Bach/ #Lambda ______ <  %.3f    cm",                        fCascSels[4] ));
                fPaveTextBookKeeping->AddText( Form("  - cos(Casc. pointing angle) ____ >  %.3f ",                               fCascSels[5] ));
                fPaveTextBookKeeping->AddText( Form("  - R_{transv}(Casc. decay) ______ >  %.3f         cm",                     fCascSels[6] ));
                fPaveTextBookKeeping->AddText( Form("  - R_{transv}(Casc. decay) ______ <  %.3f     cm",                         fCascSels[7] ));
        }
        else{   fPaveTextBookKeeping->AddText("A. No rerunning of the V0/Casc. vertexers ... See std cuts in (AliRoot+Rec.C) used for this prod. cycle");}

        fPaveTextBookKeeping->AddText("- - - - - - - - - - - ");
        
        if(fkQualityCutZprimVtxPos)      fPaveTextBookKeeping->AddText("B. Quality Cut(prim. Vtx z-Pos)    = ON  ");
        else                             fPaveTextBookKeeping->AddText("B. Quality Cut(prim. Vtx z-Pos)    = Off ");
        if(fkQualityCutNoTPConlyPrimVtx) fPaveTextBookKeeping->AddText("C. Quality Cut(No TPC prim. vtx) = ON  ");
        else                             fPaveTextBookKeeping->AddText("C. Quality Cut(No TPC prim. vtx) = Off ");
        if(fkQualityCutTPCrefit)         fPaveTextBookKeeping->AddText("D. Quality Cut(TPCrefit)               = ON  ");
        else                             fPaveTextBookKeeping->AddText("D. Quality Cut(TPCrefit)               = Off ");
        if(fkQualityCut80TPCcls)         fPaveTextBookKeeping->AddText("E. Quality Cut(80 TPC clusters)   = ON  ");
        else                             fPaveTextBookKeeping->AddText("E. Quality Cut(80 TPC clusters)   = Off ");
        if(fkExtraSelections)            fPaveTextBookKeeping->AddText("F. Extra Analysis Selections         = ON  ");
        else                             fPaveTextBookKeeping->AddText("F. Extra Analysis Selections         = Off ");

        fPaveTextBookKeeping->AddText("- - - - - - - - - - - ");

        fPaveTextBookKeeping->AddText("G. TPC Aleph Param : ");
        fPaveTextBookKeeping->AddText( Form("   - fAlephParam [0] =  %.5g", fAlephParameters[0] ));
        fPaveTextBookKeeping->AddText( Form("   - fAlephParam [1] =  %.5g", fAlephParameters[1] ));
        fPaveTextBookKeeping->AddText( Form("   - fAlephParam [2] =  %.5g", fAlephParameters[2] ));
        fPaveTextBookKeeping->AddText( Form("   - fAlephParam [3] =  %.5g", fAlephParameters[3] ));
        fPaveTextBookKeeping->AddText( Form("   - fAlephParam [4] =  %.5g", fAlephParameters[4] ));
}       
*/ 
 
	// - General histos
	//--------------
 
if(! fHistCascadeMultiplicityBeforeTrigSel) {
	if(fCollidingSystems)// AA collisions
		fHistCascadeMultiplicityBeforeTrigSel = new TH1F("fHistCascadeMultiplicityBeforeTrigSel", 
			"Cascades per event (before Trig. Sel.);Nbr of Cascades/Evt;Events", 
			100, 0, 100); 		
	else // pp collisions
	  	fHistCascadeMultiplicityBeforeTrigSel = new TH1F("fHistCascadeMultiplicityBeforeTrigSel", 
			"Cascades per event (before Trig. Sel.);Nbr of Cascades/Evt;Events", 
			25, 0, 25);
	fListHistCascade->Add(fHistCascadeMultiplicityBeforeTrigSel);
}

        // - Histos for events passing the trigger selection
        //--------------
        
if(! fHistCascadeMultiplicityForTrigEvt) {
	if(fCollidingSystems)// AA collisions
		fHistCascadeMultiplicityForTrigEvt = new TH1F("fHistCascadeMultiplicityForTrigEvt", 
			"Cascades per event (for triggered evt);Nbr of Cascades/Evt;Events", 
			100, 0, 100); 		
	else // pp collisions
	  	fHistCascadeMultiplicityForTrigEvt = new TH1F("fHistCascadeMultiplicityForTrigEvt", 
			"Cascades per event (for triggered evt);Nbr of Cascades/Evt;Events", 
			25, 0, 25);
	fListHistCascade->Add(fHistCascadeMultiplicityForTrigEvt);
}

 
if(! fHistTrackMultiplicityForTrigEvt) {	
	if(fCollidingSystems)// AA collisions	
		fHistTrackMultiplicityForTrigEvt = new TH1F("fHistTrackMultiplicityForTrigEvt", 
			"Track Multiplicity (for triggered evt);Nbr of tracks/Evt;Events", 
			200, 0, 20000); 		
	else // pp collisions
		fHistTrackMultiplicityForTrigEvt = new TH1F("fHistTrackMultiplicityForTrigEvt", 
			"Track Multiplicity (for triggered evt);Nbr of tracks/Evt;Events", 
			300, 0, 300);
	fListHistCascade->Add(fHistTrackMultiplicityForTrigEvt);
}

if(! fHistTPCrefitTrackMultiplicityForTrigEvt) {	
	if(fCollidingSystems)// AA collisions	
		fHistTPCrefitTrackMultiplicityForTrigEvt = new TH1F("fHistTPCrefitTrackMultiplicityForTrigEvt", 
			"TPCrefit track Multiplicity (for triggered evt);Nbr of TPCrefit tracks/Evt;Events", 
			200, 0, 20000); 		
	else // pp collisions
		fHistTPCrefitTrackMultiplicityForTrigEvt = new TH1F("fHistTPCrefitTrackMultiplicityForTrigEvt", 
			"TPCrefit track Multiplicity (for triggered evt);Nbr of TPCrefit tracks/Evt;Events", 
			300, 0, 300);
	fListHistCascade->Add(fHistTPCrefitTrackMultiplicityForTrigEvt);
}


if(! fHistPrimaryTrackMultiplicityForTrigEvt) {	
	if(fCollidingSystems)// AA collisions	
		fHistPrimaryTrackMultiplicityForTrigEvt = new TH1F("fHistPrimaryTrackMultiplicityForTrigEvt", 
			"Primary track Multiplicity (for triggered evt);Nbr of primary tracks/Evt;Events", 
			100, 0, 10000); 		
	else // pp collisions
		fHistPrimaryTrackMultiplicityForTrigEvt = new TH1F("fHistPrimaryTrackMultiplicityForTrigEvt", 
			"Primary track Multiplicity (for triggered evt);Nbr of primary tracks/Evt;Events", 
			200, 0, 200);
	fListHistCascade->Add(fHistPrimaryTrackMultiplicityForTrigEvt);
}


        // - Histos for events passing the trigger selection + |z(prim. vertex)| < XX cm
        //--------------
        
if(! fHistCascadeMultiplicityForTrigEvtAndZprimVtx) {
        if(fCollidingSystems)// AA collisions
                fHistCascadeMultiplicityForTrigEvtAndZprimVtx = new TH1F("fHistCascadeMultiplicityForTrigEvtAndZprimVtx", 
                        "Cascades per event;Nbr of Cascades/Evt;Events", 
                        100, 0, 100);           
        else // pp collisions
                fHistCascadeMultiplicityForTrigEvtAndZprimVtx = new TH1F("fHistCascadeMultiplicityForTrigEvtAndZprimVtx", 
                        "Cascades per event;Nbr of Cascades/Evt;Events", 
                        25, 0, 25);
        fListHistCascade->Add(fHistCascadeMultiplicityForTrigEvtAndZprimVtx);
}


        // - Histos for events passing the trigger selection + |z(prim. vertex)| < XX cm + after pile-up rejection
        //--------------
        
if(! fHistCascadeMultiplicityForTrigEvtNonPiledUpAndZprimVtx) {
        if(fCollidingSystems)// AA collisions
                fHistCascadeMultiplicityForTrigEvtNonPiledUpAndZprimVtx = new TH1F("fHistCascadeMultiplicityForTrigEvtNonPiledUpAndZprimVtx", 
                        "Cascades per event;Nbr of Cascades/Evt;Events", 
                        100, 0, 100);           
        else // pp collisions
                fHistCascadeMultiplicityForTrigEvtNonPiledUpAndZprimVtx = new TH1F("fHistCascadeMultiplicityForTrigEvtNonPiledUpAndZprimVtx", 
                        "Cascades per event;Nbr of Cascades/Evt;Events", 
                        25, 0, 25);
        fListHistCascade->Add(fHistCascadeMultiplicityForTrigEvtNonPiledUpAndZprimVtx);
}


        // - Histos for events passing the event selection at the analysis level
        //--------------
        
if(! fHistCascadeMultiplicityForSelEvt) {
	if(fCollidingSystems)// AA collisions
		fHistCascadeMultiplicityForSelEvt = new TH1F("fHistCascadeMultiplicityForSelEvt", 
			"Cascades per event;Nbr of Cascades/Evt;Events", 
			100, 0, 100); 		
	else // pp collisions
	  	fHistCascadeMultiplicityForSelEvt = new TH1F("fHistCascadeMultiplicityForSelEvt", 
			"Cascades per event;Nbr of Cascades/Evt;Events", 
			25, 0, 25);
	fListHistCascade->Add(fHistCascadeMultiplicityForSelEvt);
}


if(! fHistPosBestPrimaryVtxXForSelEvt ){
	fHistPosBestPrimaryVtxXForSelEvt   = new TH1F( "fHistPosBestPrimaryVtxXForSelEvt" , "Best Prim. Vertex Position in x; x (cm); Events" , 360, -0.9, 0.9 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxXForSelEvt);
}

if(! fHistPosBestPrimaryVtxYForSelEvt){
	fHistPosBestPrimaryVtxYForSelEvt   = new TH1F( "fHistPosBestPrimaryVtxYForSelEvt" , "Best Prim. Vertex Position in y; y (cm); Events" , 360, -0.9, 0.9 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxYForSelEvt);
}

if(! fHistPosBestPrimaryVtxZForSelEvt ){
	fHistPosBestPrimaryVtxZForSelEvt   = new TH1F( "fHistPosBestPrimaryVtxZForSelEvt" , "Best Prim. Vertex Position in z; z (cm); Events" , 300, -30.0, 30.0 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxZForSelEvt);
}




        // - Histos for events containing at least ONE CASCADE
        //--------------
        
if(! fHistTPCrefitTrackMultiplicityForCascadeEvt) {
	if(fCollidingSystems)// AA collisions	
		fHistTPCrefitTrackMultiplicityForCascadeEvt = new TH1F("fHistTPCrefitTrackMultiplicityForCascadeEvt", 
			"TPCrefit track Multiplicity (for evt with Casc.);Nbr of TPCrefit tracks/Evt with cascade(s);Events", 
			200, 0, 20000); 		
	else // pp collisions
		fHistTPCrefitTrackMultiplicityForCascadeEvt = new TH1F("fHistTPCrefitTrackMultiplicityForCascadeEvt", 
			"TPCrefit track Multiplicity (for evt with Casc.);Nbr of TPCrefit tracks/Evt with cascade(s);Events", 
			300, 0, 300);
	fListHistCascade->Add(fHistTPCrefitTrackMultiplicityForCascadeEvt);
}

if(! fHistPrimaryTrackMultiplicityForCascadeEvt) {	
	if(fCollidingSystems)// AA collisions	
		fHistPrimaryTrackMultiplicityForCascadeEvt = new TH1F("fHistPrimaryTrackMultiplicityForCascadeEvt", 
			"Primary track Multiplicity (for evt with Casc.);Nbr of primary tracks/Evt;Events", 
			100, 0, 10000); 		
	else // pp collisions
		fHistPrimaryTrackMultiplicityForCascadeEvt = new TH1F("fHistPrimaryTrackMultiplicityForCascadeEvt", 
			"Primary track Multiplicity (for evt with Casc.);Nbr of primary tracks/Evt;Events", 
			200, 0, 200);
	fListHistCascade->Add(fHistPrimaryTrackMultiplicityForCascadeEvt);
}

if(! fHistPosV0TPCClusters ){
        fHistPosV0TPCClusters = new TH1F("fHistPosV0TPCClusters", "TPC clusters for Pos. V0 daughter track, in Casc; Nbr of TPC clusters (V0 Pos.); Track counts", 165, 0.0 ,165.0);
        fListHistCascade->Add(fHistPosV0TPCClusters);
}

if(! fHistNegV0TPCClusters ){
        fHistNegV0TPCClusters = new TH1F("fHistNegV0TPCClusters", "TPC clusters for Neg. V0 daughter track, in Casc; Nbr of TPC clusters (V0 Neg.); Track counts", 165, 0.0 ,165.0);
        fListHistCascade->Add(fHistNegV0TPCClusters);
}

if(! fHistBachTPCClusters ){
        fHistBachTPCClusters = new TH1F("fHistBachTPCClusters", "TPC clusters for Bachelor track; Nbr of TPC clusters (Bach); Track counts", 165, 0.0 ,165.0);
        fListHistCascade->Add(fHistBachTPCClusters);
}







if(! fHistVtxStatus ){
	fHistVtxStatus   = new TH1F( "fHistVtxStatus" , "Does a Trckg Prim.vtx exist ?; true=1 or false=0; Nb of Events" , 4, -1.0, 3.0 );
	fListHistCascade->Add(fHistVtxStatus);
}


	// - Vertex Positions
  
if(! fHistPosTrkgPrimaryVtxXForCascadeEvt ){
	fHistPosTrkgPrimaryVtxXForCascadeEvt   = new TH1F( "fHistPosTrkgPrimaryVtxXForCascadeEvt" , "Trkg Prim. Vertex Position in x; x (cm); Events" , 360, -0.9, 0.9 );
	fListHistCascade->Add(fHistPosTrkgPrimaryVtxXForCascadeEvt);
}


if(! fHistPosTrkgPrimaryVtxYForCascadeEvt){
	fHistPosTrkgPrimaryVtxYForCascadeEvt   = new TH1F( "fHistPosTrkgPrimaryVtxYForCascadeEvt" , "Trkg Prim. Vertex Position in y; y (cm); Events" , 360, -0.9, 0.9 );
	fListHistCascade->Add(fHistPosTrkgPrimaryVtxYForCascadeEvt);
}

if(! fHistPosTrkgPrimaryVtxZForCascadeEvt ){
	fHistPosTrkgPrimaryVtxZForCascadeEvt   = new TH1F( "fHistPosTrkgPrimaryVtxZForCascadeEvt" , "Trkg Prim. Vertex Position in z; z (cm); Events" , 300, -30.0, 30.0 );
	fListHistCascade->Add(fHistPosTrkgPrimaryVtxZForCascadeEvt);
}

if(! fHistTrkgPrimaryVtxRadius ){
	fHistTrkgPrimaryVtxRadius  = new TH1F( "fHistTrkgPrimaryVtxRadius",  "Trkg Prim. Vertex radius; r (cm); Events" , 150, 0., 15.0 );
	fListHistCascade->Add(fHistTrkgPrimaryVtxRadius);
}




if(! fHistPosBestPrimaryVtxXForCascadeEvt ){
	fHistPosBestPrimaryVtxXForCascadeEvt   = new TH1F( "fHistPosBestPrimaryVtxXForCascadeEvt" , "Best Prim. Vertex Position in x; x (cm); Events" , 360, -0.9, 0.9 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxXForCascadeEvt);
}

if(! fHistPosBestPrimaryVtxYForCascadeEvt){
	fHistPosBestPrimaryVtxYForCascadeEvt   = new TH1F( "fHistPosBestPrimaryVtxYForCascadeEvt" , "Best Prim. Vertex Position in y; y (cm); Events" , 360, -0.9, 0.9 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxYForCascadeEvt);
}

if(! fHistPosBestPrimaryVtxZForCascadeEvt ){
	fHistPosBestPrimaryVtxZForCascadeEvt   = new TH1F( "fHistPosBestPrimaryVtxZForCascadeEvt" , "Best Prim. Vertex Position in z; z (cm); Events" , 300, -30.0, 30.0 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxZForCascadeEvt);
}

if(! fHistBestPrimaryVtxRadius ){
	fHistBestPrimaryVtxRadius  = new TH1F( "fHistBestPrimaryVtxRadius",  "Best Prim.  vertex radius; r (cm); Events" , 150, 0., 15.0 );
	fListHistCascade->Add(fHistBestPrimaryVtxRadius);
}

if(! f2dHistTrkgPrimVtxVsBestPrimVtx) {
	f2dHistTrkgPrimVtxVsBestPrimVtx = new TH2F( "f2dHistTrkgPrimVtxVsBestPrimVtx", "r_{Trck Prim. Vtx} Vs r_{Best Prim. Vtx}; r_{Track Vtx} (cm); r_{Best Vtx} (cm)", 300, 0., 15.0, 300, 0., 15.);
	fListHistCascade->Add(f2dHistTrkgPrimVtxVsBestPrimVtx);
}




// - Typical histos for cascades


if(! fHistEffMassXi) {
     fHistEffMassXi = new TH1F("fHistEffMassXi", "Cascade candidates ; Invariant Mass (GeV/c^{2}) ; Counts", 400, 1.2, 2.0);
     fListHistCascade->Add(fHistEffMassXi);
}
   
if(! fHistChi2Xi ){
	fHistChi2Xi = new TH1F("fHistChi2Xi", "Cascade #chi^{2}; #chi^{2}; Number of Cascades", 160, 0, 40);
	fListHistCascade->Add(fHistChi2Xi);
}
  
if(! fHistDcaXiDaughters ){
	fHistDcaXiDaughters = new TH1F( "fHistDcaXiDaughters",  "DCA between Xi Daughters; DCA (cm) ; Number of Cascades", 100, 0., 0.5);
	fListHistCascade->Add(fHistDcaXiDaughters);
}

if(! fHistDcaBachToPrimVertex) {
	fHistDcaBachToPrimVertex = new TH1F("fHistDcaBachToPrimVertex", "DCA of Bach. to Prim. Vertex;DCA (cm);Number of Cascades", 250, 0., 0.25);
	fListHistCascade->Add(fHistDcaBachToPrimVertex);
}

if(! fHistXiCosineOfPointingAngle) {
	fHistXiCosineOfPointingAngle = new TH1F("fHistXiCosineOfPointingAngle", "Cosine of Xi Pointing Angle; Cos (Xi Point.Angl);Number of Xis", 200, 0.99, 1.0);
	fListHistCascade->Add(fHistXiCosineOfPointingAngle);
}

if(! fHistXiRadius ){
	fHistXiRadius  = new TH1F( "fHistXiRadius",  "Casc. decay transv. radius; r (cm); Counts" , 1050, 0., 105.0 );
	fListHistCascade->Add(fHistXiRadius);
}


// - Histos about ~ the "V0 part" of the cascade,  coming by inheritance from AliESDv0



if (! fHistMassLambdaAsCascDghter) {
     fHistMassLambdaAsCascDghter = new TH1F("fHistMassLambdaAsCascDghter","#Lambda associated to Casc. candidates;Eff. Mass (GeV/c^{2});Counts", 300,1.00,1.3);
    fListHistCascade->Add(fHistMassLambdaAsCascDghter);
}

if (! fHistV0Chi2Xi) {
	fHistV0Chi2Xi = new TH1F("fHistV0Chi2Xi", "V0 #chi^{2}, in cascade; #chi^{2};Counts", 160, 0, 40);
	fListHistCascade->Add(fHistV0Chi2Xi);
}

if (! fHistDcaV0DaughtersXi) {
	fHistDcaV0DaughtersXi = new TH1F("fHistDcaV0DaughtersXi", "DCA between V0 daughters, in cascade;DCA (cm);Number of V0s", 120, 0., 0.6);
	fListHistCascade->Add(fHistDcaV0DaughtersXi);
}

if (! fHistDcaV0ToPrimVertexXi) {
	fHistDcaV0ToPrimVertexXi = new TH1F("fHistDcaV0ToPrimVertexXi", "DCA of V0 to Prim. Vertex, in cascade;DCA (cm);Number of Cascades", 200, 0., 1.);
	fListHistCascade->Add(fHistDcaV0ToPrimVertexXi);
}

if (! fHistV0CosineOfPointingAngleXi) {
	fHistV0CosineOfPointingAngleXi = new TH1F("fHistV0CosineOfPointingAngleXi", "Cosine of V0 Pointing Angle, in cascade;Cos(V0 Point. Angl); Counts", 200, 0.98, 1.0);
	fListHistCascade->Add(fHistV0CosineOfPointingAngleXi);
}

if (! fHistV0RadiusXi) {
	fHistV0RadiusXi = new TH1F("fHistV0RadiusXi", "V0 decay radius, in cascade; radius (cm); Counts", 1050, 0., 105.0);
	fListHistCascade->Add(fHistV0RadiusXi);
}

if (! fHistDcaPosToPrimVertexXi) {
	fHistDcaPosToPrimVertexXi = new TH1F("fHistDcaPosToPrimVertexXi", "DCA of V0 pos daughter to Prim. Vertex;DCA (cm);Counts", 300, 0, 3);
	fListHistCascade->Add(fHistDcaPosToPrimVertexXi);
}

if (! fHistDcaNegToPrimVertexXi) {
	fHistDcaNegToPrimVertexXi = new TH1F("fHistDcaNegToPrimVertexXi", "DCA of V0 neg daughter to Prim. Vertex;DCA (cm);Counts", 300, 0, 3);
	fListHistCascade->Add(fHistDcaNegToPrimVertexXi);
}




	// - Effective mass histos for cascades.
// By cascade hyp  
if (! fHistMassXiMinus) {
    fHistMassXiMinus = new TH1F("fHistMassXiMinus","#Xi^{-} candidates;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 400,1.2,2.0);
    fListHistCascade->Add(fHistMassXiMinus);
}
  
if (! fHistMassXiPlus) {
    fHistMassXiPlus = new TH1F("fHistMassXiPlus","#Xi^{+} candidates;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",400,1.2,2.0);
    fListHistCascade->Add(fHistMassXiPlus);
}

if (! fHistMassOmegaMinus) {
	fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus","#Omega^{-} candidates;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassOmegaMinus);
}
 
if (! fHistMassOmegaPlus) {
	fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus","#Omega^{+} candidates;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassOmegaPlus);
}

// By cascade hyp + bachelor PID
if (! fHistMassWithCombPIDXiMinus) {
    fHistMassWithCombPIDXiMinus = new TH1F("fHistMassWithCombPIDXiMinus","#Xi^{-} candidates, with Bach. comb. PID;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 400,1.2,2.0);
    fListHistCascade->Add(fHistMassWithCombPIDXiMinus);
}
  
if (! fHistMassWithCombPIDXiPlus) {
    fHistMassWithCombPIDXiPlus = new TH1F("fHistMassWithCombPIDXiPlus","#Xi^{+} candidates, with Bach. comb. PID;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",400,1.2,2.0);
    fListHistCascade->Add(fHistMassWithCombPIDXiPlus);
}

if (! fHistMassWithCombPIDOmegaMinus) {
	fHistMassWithCombPIDOmegaMinus = new TH1F("fHistMassWithCombPIDOmegaMinus","#Omega^{-} candidates, with Bach. comb. PID;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassWithCombPIDOmegaMinus);
}
 
if (! fHistMassWithCombPIDOmegaPlus) {
	fHistMassWithCombPIDOmegaPlus = new TH1F("fHistMassWithCombPIDOmegaPlus","#Omega^{+} candidates, with Bach. comb. PID;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassWithCombPIDOmegaPlus);
}



	// - Complements for QA

if(! fHistXiTransvMom ){
	fHistXiTransvMom  = new TH1F( "fHistXiTransvMom" , "#Xi transverse momentum (cand. around the mass peak) ; p_{t}(#Xi) (GeV/c); Counts", 100, 0.0, 10.0);
	fListHistCascade->Add(fHistXiTransvMom);
}

if(! fHistXiTotMom ){
	fHistXiTotMom  = new TH1F( "fHistXiTotMom" , "#Xi momentum norm (cand. around the mass peak); p_{tot}(#Xi) (GeV/c); Counts", 150, 0.0, 15.0);
	fListHistCascade->Add(fHistXiTotMom);
}


if(! fHistBachTransvMomXi ){
	fHistBachTransvMomXi  = new TH1F( "fHistBachTransvMomXi" , "#Xi Bach. transverse momentum (cand. around the mass peak) ; p_{t}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
	fListHistCascade->Add(fHistBachTransvMomXi);
}

if(! fHistBachTotMomXi ){
	fHistBachTotMomXi  = new TH1F( "fHistBachTotMomXi" , "#Xi Bach. momentum norm (cand. around the mass peak); p_{tot}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
	fListHistCascade->Add(fHistBachTotMomXi);
}


if(! fHistChargeXi ){
	fHistChargeXi  = new TH1F( "fHistChargeXi" , "Charge of casc. candidates ; Sign ; Counts", 5, -2.0, 3.0);
	fListHistCascade->Add(fHistChargeXi);
}


if (! fHistV0toXiCosineOfPointingAngle) {
	fHistV0toXiCosineOfPointingAngle = new TH1F("fHistV0toXiCosineOfPointingAngle", "Cos. of V0 Ptng Angl / Xi vtx ;Cos(V0 Point. Angl / Xi vtx); Counts", 100, 0.99, 1.0);
	fListHistCascade->Add(fHistV0toXiCosineOfPointingAngle);
}


if(! fHistRapXi ){
	fHistRapXi  = new TH1F( "fHistRapXi" , "Rapidity of #Xi candidates (around the mass peak); y ; Counts", 200, -5.0, 5.0);
	fListHistCascade->Add(fHistRapXi);
}

if(! fHistRapOmega ){
	fHistRapOmega  = new TH1F( "fHistRapOmega" , "Rapidity of #Omega candidates (around the mass peak); y ; Counts", 200, -5.0, 5.0);
	fListHistCascade->Add(fHistRapOmega);
}

if(! fHistEtaXi ){
	fHistEtaXi  = new TH1F( "fHistEtaXi" , "Pseudo-rap. of #Xi candidates (around the mass peak) ; #eta ; Counts", 120, -3.0, 3.0);
	fListHistCascade->Add(fHistEtaXi);
}

if(! fHistThetaXi ){
	fHistThetaXi  = new TH1F( "fHistThetaXi" , "#theta of #Xi candidates (around the mass peak); #theta (deg) ; Counts", 180, 0., 180.0);
	fListHistCascade->Add(fHistThetaXi);
}

if(! fHistPhiXi ){
	fHistPhiXi  = new TH1F( "fHistPhiXi" , "#phi of #Xi candidates (around the mass peak); #phi (deg) ; Counts", 360, 0., 360.);
	fListHistCascade->Add(fHistPhiXi);
}


if(! fHistcTauXiMinus){
        fHistcTauXiMinus = new TH1F("fHistcTauXiMinus", "Lifetime c.#tau for #Xi^{-}; L_{3D}.m_{PDG}(#Xi^{-}) / p_{3D} (cm); Counts", 100, 0., 50.);
        fListHistCascade->Add(fHistcTauXiMinus);
}

if(! fHistcTauXiPlus){
        fHistcTauXiPlus = new TH1F("fHistcTauXiPlus", "Lifetime c.#tau for #Xi^{+}; L_{3D}.m_{PDG}(#bar{#Xi}^{+}) / p_{3D} (cm); Counts", 100, 0., 50.);
        fListHistCascade->Add(fHistcTauXiPlus);
}

if(! fHistcTauOmegaMinus){
        fHistcTauOmegaMinus = new TH1F("fHistcTauOmegaMinus", "Lifetime c.#tau for #Omega^{-}; L_{3D}.m_{PDG}(#Omega^{-}) / p_{3D} (cm); Counts", 100, 0., 50.);
        fListHistCascade->Add(fHistcTauOmegaMinus);
}

if(! fHistcTauOmegaPlus){
        fHistcTauOmegaPlus = new TH1F("fHistcTauOmegaPlus", "Lifetime c.#tau for #Omega^{+}; L_{3D}.m_{PDG}(#bar{#Omega}^{+}) / p_{3D} (cm); Counts", 100, 0., 50.);
        fListHistCascade->Add(fHistcTauOmegaPlus);
}


if(! f2dHistArmenteros) {
	f2dHistArmenteros = new TH2F( "f2dHistArmenteros", "#alpha_{Arm}(casc. cand.) Vs Pt_{Arm}(casc. cand.); #alpha_{Arm} ; Pt_{Arm} (GeV/c)", 140, -1.2, 1.2, 300, 0., 0.3);
	fListHistCascade->Add(f2dHistArmenteros);
}

//-------

if(! f2dHistEffMassLambdaVsEffMassXiMinus) {
	f2dHistEffMassLambdaVsEffMassXiMinus = new TH2F( "f2dHistEffMassLambdaVsEffMassXiMinus", "M_{#Lambda} Vs M_{#Xi^{-} candidates} ; Inv. M_{#Lambda^{0}} (GeV/c^{2}) ; M( #Lambda , #pi^{-} ) (GeV/c^{2})", 300, 1.1,1.13, 400, 1.2, 2.0);
	fListHistCascade->Add(f2dHistEffMassLambdaVsEffMassXiMinus);
}

if(! f2dHistEffMassXiVsEffMassOmegaMinus) {
	f2dHistEffMassXiVsEffMassOmegaMinus = new TH2F( "f2dHistEffMassXiVsEffMassOmegaMinus", "M_{#Xi^{-} candidates} Vs M_{#Omega^{-} candidates} ; M( #Lambda , #pi^{-} ) (GeV/c^{2}) ; M( #Lambda , K^{-} ) (GeV/c^{2})", 400, 1.2, 2.0, 500, 1.5, 2.5);
	fListHistCascade->Add(f2dHistEffMassXiVsEffMassOmegaMinus);
}

if(! f2dHistEffMassLambdaVsEffMassXiPlus) {
	f2dHistEffMassLambdaVsEffMassXiPlus = new TH2F( "f2dHistEffMassLambdaVsEffMassXiPlus", "M_{#Lambda} Vs M_{#Xi^{+} candidates} ; Inv. M_{#Lambda^{0}} (GeV/c^{2}) ; M( #Lambda , #pi^{+} ) (GeV/c^{2})", 300, 1.1,1.13, 400, 1.2, 2.0);
	fListHistCascade->Add(f2dHistEffMassLambdaVsEffMassXiPlus);
}

if(! f2dHistEffMassXiVsEffMassOmegaPlus) {
	f2dHistEffMassXiVsEffMassOmegaPlus = new TH2F( "f2dHistEffMassXiVsEffMassOmegaPlus", "M_{#Xi^{+} candidates} Vs M_{#Omega^{+} candidates} ; M( #Lambda , #pi^{+} ) (GeV/c^{2}) ; M( #Lambda , K^{+} ) (GeV/c^{2})", 400, 1.2, 2.0, 500, 1.5, 2.5);
	fListHistCascade->Add(f2dHistEffMassXiVsEffMassOmegaPlus);
}

//-------

if(! f2dHistXiRadiusVsEffMassXiMinus) {
	f2dHistXiRadiusVsEffMassXiMinus = new TH2F( "f2dHistXiRadiusVsEffMassXiMinus", "Transv. R_{Xi Decay} Vs M_{#Xi^{-} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{-} ) (GeV/c^{2}) ", 450, 0., 45.0, 400, 1.2, 2.0);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiMinus);
}

if(! f2dHistXiRadiusVsEffMassXiPlus) {
	f2dHistXiRadiusVsEffMassXiPlus = new TH2F( "f2dHistXiRadiusVsEffMassXiPlus", "Transv. R_{Xi Decay} Vs M_{#Xi^{+} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{+} ) (GeV/c^{2}) ", 450, 0., 45.0, 400, 1.2, 2.0);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiPlus);
}

if(! f2dHistXiRadiusVsEffMassOmegaMinus) {
	f2dHistXiRadiusVsEffMassOmegaMinus = new TH2F( "f2dHistXiRadiusVsEffMassOmegaMinus", "Transv. R_{Xi Decay} Vs M_{#Omega^{-} candidates}; r_{cascade} (cm); M( #Lambda , K^{-} ) (GeV/c^{2}) ", 450, 0., 45.0, 500, 1.5, 2.5);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassOmegaMinus);
}

if(! f2dHistXiRadiusVsEffMassOmegaPlus) {
	f2dHistXiRadiusVsEffMassOmegaPlus = new TH2F( "f2dHistXiRadiusVsEffMassOmegaPlus", "Transv. R_{Xi Decay} Vs M_{#Omega^{+} candidates}; r_{cascade} (cm); M( #Lambda , K^{+} ) (GeV/c^{2}) ", 450, 0., 45.0, 500, 1.5, 2.5);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassOmegaPlus);
}

//------

if(! f2dHistTPCdEdxOfCascDghters){
        f2dHistTPCdEdxOfCascDghters = new TH2F( "f2dHistTPCdEdxOfCascDghters", "TPC dE/dx of the cascade daughters; charge x || #vec{p}_{TPC inner wall}(Casc. daughter) || (GeV/c); TPC signal (ADC) ", 2000, -10.0, 10.0, 450, 0., 900.);
	fListHistCascade->Add(f2dHistTPCdEdxOfCascDghters);
}



// Part 2 : Raw material for yield extraction -------

if(! f3dHistXiPtVsEffMassVsYXiMinus) {
	f3dHistXiPtVsEffMassVsYXiMinus = new TH3F( "f3dHistXiPtVsEffMassVsYXiMinus", "Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}; Pt_{cascade} (GeV/c); M( #Lambda , #pi^{-} ) (GeV/c^{2}) ;Y_{#Xi} ", 100, 0., 10.0, 400, 1.2, 2.0, 44, -1.1,1.1);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYXiMinus);
}

if(! f3dHistXiPtVsEffMassVsYXiPlus) {
	f3dHistXiPtVsEffMassVsYXiPlus = new TH3F( "f3dHistXiPtVsEffMassVsYXiPlus", "Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}; Pt_{cascade} (GeV/c); M( #Lambda , #pi^{+} ) (GeV/c^{2}); Y_{#Xi}", 100, 0., 10.0, 400, 1.2, 2.0, 44, -1.1,1.1);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYXiPlus);
}

if(! f3dHistXiPtVsEffMassVsYOmegaMinus) {
	f3dHistXiPtVsEffMassVsYOmegaMinus = new TH3F( "f3dHistXiPtVsEffMassVsYOmegaMinus", "Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{-} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 44, -1.1,1.1);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYOmegaMinus);
}

if(! f3dHistXiPtVsEffMassVsYOmegaPlus) {
	f3dHistXiPtVsEffMassVsYOmegaPlus = new TH3F( "f3dHistXiPtVsEffMassVsYOmegaPlus", "Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{+} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 44, -1.1,1.1);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYOmegaPlus);
}

//--
if(! fCFContCascadePIDXiMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 400;
  lNbBinsPerVar[2] = 44;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDXiMinus = new AliCFContainer("fCFContCascadePIDXiMinus","Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDXiMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDXiMinus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
  fCFContCascadePIDXiMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  if(fCollidingSystems) 
        fCFContCascadePIDXiMinus->SetBinLimits(3, 0.0, 20000.0  );    // Primary track Multiplicity
  else
        fCFContCascadePIDXiMinus->SetBinLimits(3, 0.0, 250.0  );     // Primary track Multiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDXiMinus->SetStepTitle(0, "No PID");
  fCFContCascadePIDXiMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDXiMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDXiMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDXiMinus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDXiMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDXiMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDXiMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDXiMinus->SetVarTitle(1, "M( #Lambda , #pi^{-} ) (GeV/c^{2})");
  fCFContCascadePIDXiMinus->SetVarTitle(2, "Y_{#Xi}");
  fCFContCascadePIDXiMinus->SetVarTitle(3, "Primary track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDXiMinus);
  
}

if(! fCFContCascadePIDXiPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 400;
  lNbBinsPerVar[2] = 44;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDXiPlus = new AliCFContainer("fCFContCascadePIDXiPlus","Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDXiPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDXiPlus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
  fCFContCascadePIDXiPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  if(fCollidingSystems) 
        fCFContCascadePIDXiPlus->SetBinLimits(3, 0.0, 20000.0  );    // Primary track Multiplicity
  else
        fCFContCascadePIDXiPlus->SetBinLimits(3, 0.0, 250.0  );     // Primary track Multiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDXiPlus->SetStepTitle(0, "No PID");
  fCFContCascadePIDXiPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDXiPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDXiPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDXiPlus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDXiPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDXiPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDXiPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDXiPlus->SetVarTitle(1, "M( #Lambda , #pi^{+} ) (GeV/c^{2})");
  fCFContCascadePIDXiPlus->SetVarTitle(2, "Y_{#Xi}");
  fCFContCascadePIDXiPlus->SetVarTitle(3, "Primary track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDXiPlus);
  
}


if(! fCFContCascadePIDOmegaMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 500;
  lNbBinsPerVar[2] = 44;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDOmegaMinus = new AliCFContainer("fCFContCascadePIDOmegaMinus","Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDOmegaMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDOmegaMinus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
  fCFContCascadePIDOmegaMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  if(fCollidingSystems) 
        fCFContCascadePIDOmegaMinus->SetBinLimits(3, 0.0, 20000.0  );    //Primary track Multiplicity
  else
        fCFContCascadePIDOmegaMinus->SetBinLimits(3, 0.0, 250.0  );     // Primary track Multiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDOmegaMinus->SetStepTitle(0, "No PID");
  fCFContCascadePIDOmegaMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDOmegaMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDOmegaMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDOmegaMinus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDOmegaMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDOmegaMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDOmegaMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDOmegaMinus->SetVarTitle(1, "M( #Lambda , K^{-} ) (GeV/c^{2})");
  fCFContCascadePIDOmegaMinus->SetVarTitle(2, "Y_{#Omega}");
  fCFContCascadePIDOmegaMinus->SetVarTitle(3, "Primary track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDOmegaMinus);
  
}

if(! fCFContCascadePIDOmegaPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 500;
  lNbBinsPerVar[2] = 44;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDOmegaPlus = new AliCFContainer("fCFContCascadePIDOmegaPlus","Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDOmegaPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDOmegaPlus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
  fCFContCascadePIDOmegaPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  if(fCollidingSystems) 
        fCFContCascadePIDOmegaPlus->SetBinLimits(3, 0.0, 20000.0  );    // Primary track Multiplicity
  else
        fCFContCascadePIDOmegaPlus->SetBinLimits(3, 0.0, 250.0  );     // Primary track Multiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDOmegaPlus->SetStepTitle(0, "No PID");
  fCFContCascadePIDOmegaPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDOmegaPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDOmegaPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDOmegaPlus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDOmegaPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDOmegaPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDOmegaPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDOmegaPlus->SetVarTitle(1, "M( #Lambda , K^{+} ) (GeV/c^{2})");
  fCFContCascadePIDOmegaPlus->SetVarTitle(2, "Y_{#Omega}");
  fCFContCascadePIDOmegaPlus->SetVarTitle(3, "Primary track Multiplicity"); 
  
  fListHistCascade->Add(fCFContCascadePIDOmegaPlus);
  
}





// Part 3 : Towards the optimisation of topological selections -------
if(! fCFContCascadeCuts){
   
	// Container meant to store all the relevant distributions corresponding to the cut variables.
	// So far, 20 variables have been identified.
	// The following will be done in quite a brut force way ... 
	// FIXME Improvement expected later (before Pb-Pb data at least)
        //          - Define a user binning to have less bins in each dimension
        //          - boolean for enabling/disbaling this CFContainer
  const	Int_t  lNbSteps      =  4 ;
  const Int_t  lNbVariables  =  20 ;
  
  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[20] = {0};
  lNbBinsPerVar[0]  = 25;
  lNbBinsPerVar[1]  = 25;
  lNbBinsPerVar[2]  = 20;
  lNbBinsPerVar[3]  = 40;
  lNbBinsPerVar[4]  = 30;
  lNbBinsPerVar[5]  = 25;
  
  lNbBinsPerVar[6]  = 20;
  lNbBinsPerVar[7]  = 40;
  lNbBinsPerVar[8]  = 40;
  lNbBinsPerVar[9]  = 25;
  lNbBinsPerVar[10] = 25;
  
  lNbBinsPerVar[11] = 75;  // 2-MeV/c2 bins
  lNbBinsPerVar[12] = 60;  // 2-MeV/c2 bins
  
  lNbBinsPerVar[13] = 100;
  
  lNbBinsPerVar[14] = 44; // 0.05 in rapidity units
  lNbBinsPerVar[15] = 44; // 0.05 in rapidity units
  
  lNbBinsPerVar[16] = 20;
 
  lNbBinsPerVar[17] = 50;
  lNbBinsPerVar[18] = 100;
  lNbBinsPerVar[19] = 24;
    
 fCFContCascadeCuts = new AliCFContainer("fCFContCascadeCuts","Container for Cascade cuts", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits 
   
  //0
  Double_t *lBinLim0  = new Double_t[ lNbBinsPerVar[0]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[0];i++)  lBinLim0[i]  = (Double_t)0.0   + (4.8  - 0.0 )/(lNbBinsPerVar[0]-1)  * (Double_t)i ;
        lBinLim0[ lNbBinsPerVar[0]  ] = 20.0;
  fCFContCascadeCuts -> SetBinLimits(0,  lBinLim0 );            // DcaXiDaughters : 0.0 to 5.0	
  delete [] lBinLim0;  
  //1
  Double_t *lBinLim1  = new Double_t[ lNbBinsPerVar[1]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[1];i++)   lBinLim1[i]  = (Double_t)0.0   + (0.24  - 0.0 )/(lNbBinsPerVar[1]-1)  * (Double_t)i ;
        lBinLim1[ lNbBinsPerVar[1]  ] = 100.0;
  fCFContCascadeCuts -> SetBinLimits(1,  lBinLim1 );            // DcaBachToPrimVertexXi : 0.0 to 0.25	
  delete [] lBinLim1;  
  //2 
  Double_t *lBinLim2  = new Double_t[ lNbBinsPerVar[2]+1 ];
        for(Int_t i=1; i< lNbBinsPerVar[2]+1;i++)   lBinLim2[i]  = (Double_t)0.81   + (1.0  - 0.81 )/(lNbBinsPerVar[2]-1)  * (Double_t) (i-1) ;   
        lBinLim2[0] = 0.0;
  fCFContCascadeCuts -> SetBinLimits(2,  lBinLim2 );            // XiCosineOfPointingAngle : 0.80 to 1.0	
  delete [] lBinLim2;  
  //3
  Double_t *lBinLim3  = new Double_t[ lNbBinsPerVar[3]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[3];i++)   lBinLim3[i]  = (Double_t)0.0   + (3.9  - 0.0 )/(lNbBinsPerVar[3]-1)  * (Double_t)i ;
        lBinLim3[ lNbBinsPerVar[3]  ] = 110.0;
  fCFContCascadeCuts -> SetBinLimits(3,  lBinLim3 );            // XiRadius : 0.0 to 4.0	
  delete [] lBinLim3;  
  //4
  fCFContCascadeCuts->SetBinLimits(4,    1.1  ,  1.13 );        // InvMassLambdaAsCascDghter
  //5
  Double_t *lBinLim5  = new Double_t[ lNbBinsPerVar[5]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[5];i++)   lBinLim5[i]  = (Double_t)0.0   + (4.8  - 0.0 )/(lNbBinsPerVar[5]-1)  * (Double_t)i ;
        lBinLim5[ lNbBinsPerVar[5]  ] = 20.0;
  fCFContCascadeCuts -> SetBinLimits(5,  lBinLim5 );            // DcaV0DaughtersXi : 0.0 to 5.0	
  delete [] lBinLim5;  
  
  
  //6 
  Double_t *lBinLim6  = new Double_t[ lNbBinsPerVar[6]+1 ];
        for(Int_t i=1; i< lNbBinsPerVar[6]+1 ;i++)   lBinLim6[i]  = (Double_t)0.81   + (1.0  - 0.81 )/(lNbBinsPerVar[6]-1)  * (Double_t) (i-1) ;   
        lBinLim6[0] = 0.0;
  fCFContCascadeCuts -> SetBinLimits(6,  lBinLim6 );            // V0CosineOfPointingAngleXi : 0.80 to 1.0	
  delete [] lBinLim6;  
  //7
  Double_t *lBinLim7  = new Double_t[ lNbBinsPerVar[7]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[7];i++)   lBinLim7[i]  = (Double_t)0.0   + (7.8  - 0.0 )/(lNbBinsPerVar[7]-1)  * (Double_t)i ;
        lBinLim7[ lNbBinsPerVar[7]  ] = 100.0;
  fCFContCascadeCuts -> SetBinLimits(7,  lBinLim7 );            // V0RadiusXi : 0.0 to 8.0	
  delete [] lBinLim7;  
  //8
  Double_t *lBinLim8  = new Double_t[ lNbBinsPerVar[8]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[8];i++)   lBinLim8[i]  = (Double_t)0.0   + (0.39  - 0.0 )/(lNbBinsPerVar[8]-1)  * (Double_t)i ;
        lBinLim8[ lNbBinsPerVar[8]  ] = 100.0;
  fCFContCascadeCuts -> SetBinLimits(8,  lBinLim8 );            // DcaV0ToPrimVertexXi : 0.0 to 0.40	
  delete [] lBinLim8;  
  //9
  Double_t *lBinLim9  = new Double_t[ lNbBinsPerVar[9]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[9];i++)   lBinLim9[i]  = (Double_t)0.0   + (0.24  - 0.0 )/(lNbBinsPerVar[9]-1)  * (Double_t)i ;
        lBinLim9[ lNbBinsPerVar[9]  ] = 100.0;
  fCFContCascadeCuts -> SetBinLimits(9,  lBinLim9 );            // DcaPosToPrimVertexXi : 0.0 to 0.25	
  delete [] lBinLim9; 
  //10
  Double_t *lBinLim10  = new Double_t[ lNbBinsPerVar[10]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[10];i++)   lBinLim10[i]  = (Double_t)0.0   + (0.24  - 0.0 )/(lNbBinsPerVar[10]-1)  * (Double_t)i ;
        lBinLim10[ lNbBinsPerVar[10]  ] = 100.0;
  fCFContCascadeCuts -> SetBinLimits(10,  lBinLim10 );            // DcaPosToPrimVertexXi : 0.0 to 0.25	
  delete [] lBinLim10; 
  
  //11
  fCFContCascadeCuts->SetBinLimits(11,   1.25 ,  1.40 );	// InvMassXi
  fCFContCascadeCuts->SetBinLimits(12,   1.62  , 1.74 );	// InvMassOmega
  fCFContCascadeCuts->SetBinLimits(13,   0.0  , 10.0  );	// XiTransvMom 
  fCFContCascadeCuts->SetBinLimits(14,  -1.1  ,  1.1  );	// Y(Xi)
  fCFContCascadeCuts->SetBinLimits(15,  -1.1  ,  1.1  );	// Y(Omega)
  fCFContCascadeCuts->SetBinLimits(16, -10.0  , 10.0  );	// BestPrimaryVtxPosZ
  if(fCollidingSystems){
  	fCFContCascadeCuts->SetBinLimits(17,   0.0, 10000.0  );    // nTrackPrimaryMultiplicity
  	fCFContCascadeCuts->SetBinLimits(18,   0.0, 10000.0  );    // SPDTrackletsMultiplicity
  }
  else{  
  	fCFContCascadeCuts->SetBinLimits(17,   0.0, 250.0  );     // nTrackPrimaryMultiplicity
  	fCFContCascadeCuts->SetBinLimits(18,   0.0, 200.0  );     // SPDTrackletsMultiplicity
  }
  fCFContCascadeCuts->SetBinLimits(19,  68.0  ,164.0  );	// BachTPCClusters
  
  
  // Regular binning definition (valid for v4-18-10-AN on)
  /*
  fCFContCascadeCuts->SetBinLimits(0,    0.0  ,  2.5  );        // DcaXiDaughters
  fCFContCascadeCuts->SetBinLimits(1,    0.0  ,  0.25 );	// DcaBachToPrimVertexXi
  fCFContCascadeCuts->SetBinLimits(2,    0.99 ,  1.0  );	// XiCosineOfPointingAngle
  fCFContCascadeCuts->SetBinLimits(3,    0.0  ,  4.0  );	// XiRadius
  fCFContCascadeCuts->SetBinLimits(4,    1.1  ,  1.15 );	// InvMassLambdaAsCascDghter
  fCFContCascadeCuts->SetBinLimits(5,    0.0  ,  1.0  );	// DcaV0DaughtersXi
  fCFContCascadeCuts->SetBinLimits(6,    0.98 ,  1.0  );	// V0CosineOfPointingAngleXi
  fCFContCascadeCuts->SetBinLimits(7,    0.0  , 20.0  );	// V0RadiusXi
  fCFContCascadeCuts->SetBinLimits(8,    0.0  ,  1.0  );	// DcaV0ToPrimVertexXi
  fCFContCascadeCuts->SetBinLimits(9,    0.0  ,  0.25 );	// DcaPosToPrimVertexXi
  fCFContCascadeCuts->SetBinLimits(10,   0.0  ,  0.25 );	// DcaNegToPrimVertexXi
  fCFContCascadeCuts->SetBinLimits(11,   1.25 ,  1.40 );	// InvMassXi
  fCFContCascadeCuts->SetBinLimits(12,   1.62  , 1.74 );	// InvMassOmega
  fCFContCascadeCuts->SetBinLimits(13,   0.0  , 10.0  );	// XiTransvMom 
  fCFContCascadeCuts->SetBinLimits(14,  -1.1  ,  1.1  );	// Y(Xi)
  fCFContCascadeCuts->SetBinLimits(15,  -1.1  ,  1.1  );	// Y(Omega)
  fCFContCascadeCuts->SetBinLimits(16, -10.0  , 10.0  );	// BestPrimaryVtxPosZ
  if(fCollidingSystems){
  	fCFContCascadeCuts->SetBinLimits(17,   0.0, 10000.0  );    // nTrackPrimaryMultiplicity
  	fCFContCascadeCuts->SetBinLimits(18,   0.0, 10000.0  );    // SPDTrackletsMultiplicity
  }
  else{  
  	fCFContCascadeCuts->SetBinLimits(17,   0.0, 250.0  );     // nTrackPrimaryMultiplicity
  	fCFContCascadeCuts->SetBinLimits(18,   0.0, 200.0  );     // SPDTrackletsMultiplicity
  }
  fCFContCascadeCuts->SetBinLimits(19,  25.0  ,165.0  );	// BachTPCClusters
  */
  
  
  // Setting the number of steps : one for each cascade species (Xi-, Xi+ and Omega-, Omega+)
  fCFContCascadeCuts->SetStepTitle(0, "#Xi^{-} candidates");
  fCFContCascadeCuts->SetStepTitle(1, "#bar{#Xi}^{+} candidates");
  fCFContCascadeCuts->SetStepTitle(2, "#Omega^{-} candidates");
  fCFContCascadeCuts->SetStepTitle(3, "#bar{#Omega}^{+} candidates");
  
  // Setting the variable title, per axis
  // fCFContCascadeCuts->SetVarTitle(40,  "Chi2Xi");
  fCFContCascadeCuts->SetVarTitle(0,  "Dca(XiDaughters) (cm)");
  fCFContCascadeCuts->SetVarTitle(1,  "Dca(Bach/PrimVertex) (cm)");
  fCFContCascadeCuts->SetVarTitle(2,  "cos(Xi pointing angle)");
  fCFContCascadeCuts->SetVarTitle(3,  "R_{2d}(Xi decay) (cm)");
  fCFContCascadeCuts->SetVarTitle(4,  "M_{#Lambda}(As Casc Dghter) (GeV/c^{2})");
   // fCFContCascadeCuts->SetVarTitle(40,  "V0Chi2Xi");
  fCFContCascadeCuts->SetVarTitle(5,  "Dca(V0 Daughters) in Xi (cm)");
  
  fCFContCascadeCuts->SetVarTitle(6,  "cos(V0 pointing Angle) in Casc");
  fCFContCascadeCuts->SetVarTitle(7,  "R_{2d}(V0 decay) (cm)");
  fCFContCascadeCuts->SetVarTitle(8,  "Dca(V0/PrimVertex) (cm)");
  fCFContCascadeCuts->SetVarTitle(9,  "Dca(Pos/PrimVertex) (cm)");
  fCFContCascadeCuts->SetVarTitle(10, "Dca(Neg/PrimVertex) (cm)");
  
  fCFContCascadeCuts->SetVarTitle(11, "Inv. Mass(Xi) (GeV/c^{2})");
  fCFContCascadeCuts->SetVarTitle(12, "Inv. Mass(Omega) (GeV/c^{2})");
  
  fCFContCascadeCuts->SetVarTitle(13, "pt(Casc.) (GeV/c)");
  //fCFContCascadeCuts->SetVarTitle(40, "V0toXiCosineOfPointingAngle");
  
  fCFContCascadeCuts->SetVarTitle(14, "Y(Xi)");
  fCFContCascadeCuts->SetVarTitle(15, "Y(Omega)");
  
  fCFContCascadeCuts->SetVarTitle(16, "Z-position(BestPrimVtx) (cm)");
  
  fCFContCascadeCuts->SetVarTitle(17, "Primary Track Multiplicity");
  fCFContCascadeCuts->SetVarTitle(18, "SPD tracklets Multiplicity");
  fCFContCascadeCuts->SetVarTitle(19, "Bach.TPC Clusters");
  
  fListHistCascade->Add(fCFContCascadeCuts);
}



// Part 4 : Angular correlation study -------

if(! fHnSpAngularCorrXiMinus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
	
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.30};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.34};

    TString strHnSparseTitle("");
    TString strAxisTitle[5];
    if(fAngularCorrelationType == "TrigAnyCasc-AssoAnyPrim" ){
        strHnSparseTitle = "Angular Correlation for #Xi^{-}: Trig = Casc. / Asso = all prim. tracks";
        strAxisTitle[0]  = " #Delta#phi(Casc,Track) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Casc,Track)";
        strAxisTitle[2]  = " Pt_{Casc} (GeV/c)";
        strAxisTitle[3]  = " Pt_{asso. track} (GeV/c)";
    }
    else if(fAngularCorrelationType == "TrigCascLeading-AssoAnyPrim"){
        strHnSparseTitle = "Angular Correlation for #Xi^{-}: Trig = Casc. (leading part.) / Asso = all prim. tracks";
        strAxisTitle[0]  = " #Delta#phi(Casc_{LEADING},Track) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Casc_{LEADING},Track)";
        strAxisTitle[2]  = " Pt(Casc_{LEADING}) (GeV/c)";
        strAxisTitle[3]  = " Pt_{asso. track} (GeV/c)";
    }
    else if(fAngularCorrelationType == "TrigLeadingTrck-AssoCasc"){
        strHnSparseTitle = "Angular Correlation for #Xi^{-}: Trig = leading track / Asso = any cascade";
        strAxisTitle[0]  = " #Delta#phi(Leading Track,Casc) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Leading Track,Casc)";
        strAxisTitle[2]  = " Pt(asso. Casc) (GeV/c)";
        strAxisTitle[3]  = " Pt_{Leading track} (GeV/c)";

    }
        strAxisTitle[4]  = " Eff. Inv Mass (GeV/c^{2})";

   fHnSpAngularCorrXiMinus = new THnSparseF("fHnSpAngularCorrXiMinus", strHnSparseTitle.Data(), 5, bins, xmin, xmax);
        fHnSpAngularCorrXiMinus->GetAxis(0)->SetTitle( strAxisTitle[0].Data() );
        fHnSpAngularCorrXiMinus->GetAxis(1)->SetTitle( strAxisTitle[1].Data() );
        fHnSpAngularCorrXiMinus->GetAxis(2)->SetTitle( strAxisTitle[2].Data() );
        fHnSpAngularCorrXiMinus->GetAxis(3)->SetTitle( strAxisTitle[3].Data() );
        fHnSpAngularCorrXiMinus->GetAxis(4)->SetTitle( strAxisTitle[4].Data() );
        fHnSpAngularCorrXiMinus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrXiMinus);
}

if(! fHnSpAngularCorrXiPlus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.30};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.34};
   
    TString strHnSparseTitle("");
    TString strAxisTitle[5];
    if(fAngularCorrelationType == "TrigAnyCasc-AssoAnyPrim" ){
        strHnSparseTitle = "Angular Correlation for #bar{#Xi}^{+}: Trig = Casc. / Asso = all prim. tracks";
        strAxisTitle[0]  = " #Delta#phi(Casc,Track) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Casc,Track)";
        strAxisTitle[2]  = " Pt_{Casc} (GeV/c)";
        strAxisTitle[3]  = " Pt_{asso. track} (GeV/c)";
    }
    else if(fAngularCorrelationType == "TrigCascLeading-AssoAnyPrim"){
        strHnSparseTitle = "Angular Correlation for #bar{#Xi}^{+}: Trig = Casc. (leading part.) / Asso = all prim. tracks";
        strAxisTitle[0]  = " #Delta#phi(Casc_{LEADING},Track) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Casc_{LEADING},Track)";
        strAxisTitle[2]  = " Pt(Casc_{LEADING}) (GeV/c)";
        strAxisTitle[3]  = " Pt_{asso. track} (GeV/c)";
    }
    else if(fAngularCorrelationType == "TrigLeadingTrck-AssoCasc"){
        strHnSparseTitle = "Angular Correlation for #bar{#Xi}^{+}: Trig = leading track / Asso = any cascade";
        strAxisTitle[0]  = " #Delta#phi(Leading Track,Casc) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Leading Track,Casc)";
        strAxisTitle[2]  = " Pt(asso. Casc) (GeV/c)";
        strAxisTitle[3]  = " Pt_{Leading track} (GeV/c)";

    }
        strAxisTitle[4]  = " Eff. Inv Mass (GeV/c^{2})";

   fHnSpAngularCorrXiPlus = new THnSparseF("fHnSpAngularCorrXiPlus", strHnSparseTitle.Data(), 5, bins, xmin, xmax);
	fHnSpAngularCorrXiPlus->GetAxis(0)->SetTitle( strAxisTitle[0].Data() );
	fHnSpAngularCorrXiPlus->GetAxis(1)->SetTitle( strAxisTitle[1].Data() );
	fHnSpAngularCorrXiPlus->GetAxis(2)->SetTitle( strAxisTitle[2].Data() );
	fHnSpAngularCorrXiPlus->GetAxis(3)->SetTitle( strAxisTitle[3].Data() );
	fHnSpAngularCorrXiPlus->GetAxis(4)->SetTitle( strAxisTitle[4].Data() );
	fHnSpAngularCorrXiPlus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrXiPlus);
}

if(! fHnSpAngularCorrOmegaMinus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
	
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.65};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.69};
   
    TString strHnSparseTitle("");
    TString strAxisTitle[5];
    if(fAngularCorrelationType == "TrigAnyCasc-AssoAnyPrim" ){
        strHnSparseTitle = "Angular Correlation for #Omega^{-}: Trig = Casc. / Asso = all prim. tracks";
        strAxisTitle[0]  = " #Delta#phi(Casc,Track) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Casc,Track)";
        strAxisTitle[2]  = " Pt_{Casc} (GeV/c)";
        strAxisTitle[3]  = " Pt_{asso. track} (GeV/c)";
    }
    else if(fAngularCorrelationType == "TrigCascLeading-AssoAnyPrim"){
        strHnSparseTitle = "Angular Correlation for #Omega^{-}: Trig = Casc. (leading part.) / Asso = all prim. tracks";
        strAxisTitle[0]  = " #Delta#phi(Casc_{LEADING},Track) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Casc_{LEADING},Track)";
        strAxisTitle[2]  = " Pt(Casc_{LEADING}) (GeV/c)";
        strAxisTitle[3]  = " Pt_{asso. track} (GeV/c)";
    }
    else if(fAngularCorrelationType == "TrigLeadingTrck-AssoCasc"){
        strHnSparseTitle = "Angular Correlation for #Omega^{-}: Trig = leading track / Asso = any cascade";
        strAxisTitle[0]  = " #Delta#phi(Leading Track,Casc) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Leading Track,Casc)";
        strAxisTitle[2]  = " Pt(asso. Casc) (GeV/c)";
        strAxisTitle[3]  = " Pt_{Leading track} (GeV/c)";

    }
        strAxisTitle[4]  = " Eff. Inv Mass (GeV/c^{2})";
    
   fHnSpAngularCorrOmegaMinus = new THnSparseF("fHnSpAngularCorrOmegaMinus", strHnSparseTitle.Data(), 5, bins, xmin, xmax);
	fHnSpAngularCorrOmegaMinus->GetAxis(0)->SetTitle( strAxisTitle[0].Data() );
	fHnSpAngularCorrOmegaMinus->GetAxis(1)->SetTitle( strAxisTitle[1].Data() );
	fHnSpAngularCorrOmegaMinus->GetAxis(2)->SetTitle( strAxisTitle[2].Data() );
	fHnSpAngularCorrOmegaMinus->GetAxis(3)->SetTitle( strAxisTitle[3].Data() );
	fHnSpAngularCorrOmegaMinus->GetAxis(4)->SetTitle( strAxisTitle[4].Data() );
	fHnSpAngularCorrOmegaMinus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrOmegaMinus);
}

if(! fHnSpAngularCorrOmegaPlus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.65};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.69};
   
    TString strHnSparseTitle("");
    TString strAxisTitle[5];
        if(fAngularCorrelationType == "TrigAnyCasc-AssoAnyPrim" ){
        strHnSparseTitle = "Angular Correlation for #bar{#Omega}^{+}: Trig = Casc. / Asso = all prim. tracks";
        strAxisTitle[0]  = " #Delta#phi(Casc,Track) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Casc,Track)";
        strAxisTitle[2]  = " Pt_{Casc} (GeV/c)";
        strAxisTitle[3]  = " Pt_{asso. track} (GeV/c)";
    }
    else if(fAngularCorrelationType == "TrigCascLeading-AssoAnyPrim"){
        strHnSparseTitle = "Angular Correlation for #bar{#Omega}^{+}: Trig = Casc. (leading part.) / Asso = all prim. tracks";
        strAxisTitle[0]  = " #Delta#phi(Casc_{LEADING},Track) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Casc_{LEADING},Track)";
        strAxisTitle[2]  = " Pt(Casc_{LEADING}) (GeV/c)";
        strAxisTitle[3]  = " Pt_{asso. track} (GeV/c)";
    }
    else if(fAngularCorrelationType == "TrigLeadingTrck-AssoCasc"){
        strHnSparseTitle = "Angular Correlation for #bar{#Omega}^{+}: Trig = leading track / Asso = any cascade";
        strAxisTitle[0]  = " #Delta#phi(Leading Track,Casc) (deg)";
        strAxisTitle[1]  = " #Delta#eta(Leading Track,Casc)";
        strAxisTitle[2]  = " Pt(asso. Casc) (GeV/c)";
        strAxisTitle[3]  = " Pt_{Leading track} (GeV/c)";

    }
        strAxisTitle[4]  = " Eff. Inv Mass (GeV/c^{2})";
   
   fHnSpAngularCorrOmegaPlus = new THnSparseF("fHnSpAngularCorrOmegaPlus", strHnSparseTitle.Data(), 5, bins, xmin, xmax);
	fHnSpAngularCorrOmegaPlus->GetAxis(0)->SetTitle( strAxisTitle[0].Data() );
	fHnSpAngularCorrOmegaPlus->GetAxis(1)->SetTitle( strAxisTitle[1].Data() );
	fHnSpAngularCorrOmegaPlus->GetAxis(2)->SetTitle( strAxisTitle[2].Data() );
	fHnSpAngularCorrOmegaPlus->GetAxis(3)->SetTitle( strAxisTitle[3].Data() );
	fHnSpAngularCorrOmegaPlus->GetAxis(4)->SetTitle( strAxisTitle[4].Data() );
	fHnSpAngularCorrOmegaPlus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrOmegaPlus);
}


PostData(1, fListHistCascade);
/* PostData(2, fPaveTextBookKeeping);*/
}// end UserCreateOutputObjects






//________________________________________________________________________
void AliAnalysisTaskCheckCascade::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

        AliESDEvent *lESDevent = 0x0;
        AliAODEvent *lAODevent = 0x0;
        Int_t    ncascades                      = -1;
        Int_t    nTrackMultiplicity             = -1;
        Int_t    nTrackWithTPCrefitMultiplicity = -1;
        Int_t    nTrackPrimaryMultiplicity      = -1;

        Short_t  lStatusTrackingPrimVtx         = -2;
        Double_t lTrkgPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
        Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
        Double_t lMagneticField                 = -10.;

	
	
  // Connect to the InputEvent	
  // After these lines, we should have an ESD/AOD event + the number of cascades in it.
		
  if(fAnalysisType == "ESD"){
	lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
	if (!lESDevent) {
		AliWarning("ERROR: lESDevent not available \n");
		return;
	}
        
        fHistCascadeMultiplicityBeforeTrigSel->Fill ( lESDevent->GetNumberOfCascades() );
        
        //-------------------------------------------------
        // 0 - Trigger managment
	// NOTE : Check the availability of the proper trigger 
	
        // 1st option
                //AliMCEventHandler *lmcEvtHandler  = dynamic_cast<AliMCEventHandler*>( (AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler() );
                //if( !lmcEvtHandler ){  // !0x0 = real data or !1 = there is an MC handler available (useMC = kTRUE in AnalysisTrainNew), so = data from MC
                //             if ( !( lESDevent->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL")) ) return;
                //}

        // 2nd option - Presuppose the presence of AliPhysicsSelectionTask FIXME
        
        UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
        Bool_t isSelected = 0;
        if(     fTriggerMaskType == "kMB")           isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
        else if(fTriggerMaskType == "kHighMult")     isSelected = (maskIsSelected & AliVEvent::kHighMult) == AliVEvent::kHighMult;
        else                                         isSelected = 1; // default = select anyway (use case = run without Phys Selection task)
        
        if ( ! isSelected ) { 
                PostData(1, fListHistCascade); 
                return;
        }
        
        //else Printf("Event selected ... \n");
	
        
        
        //-------------------------------------------------
        // 1 - Cascade vertexer (ESD)
        
        if(fkRerunV0CascVertexers){ // FIXME : relaunch V0 and Cascade vertexers
                if(fAnalysisType == "ESD" ){
//                         lESDevent->ResetCascades();
//                         lESDevent->ResetV0s();
//                 
//                         AliV0vertexer lV0vtxer;
//                         AliCascadeVertexer lCascVtxer;
//                 
//                         lV0vtxer.SetDefaultCuts(fV0Sels);
//                         lCascVtxer.SetDefaultCuts(fCascSels);
//                 
//                         lV0vtxer.Tracks2V0vertices(lESDevent);
//                         lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
               }
        }// end if(RelaunchV0CascVertexers)
        
        //-------------------------------------------------
	ncascades                      = lESDevent->GetNumberOfCascades();
        nTrackWithTPCrefitMultiplicity = DoESDTrackWithTPCrefitMultiplicity(lESDevent);
        nTrackPrimaryMultiplicity      = fESDtrackCuts->CountAcceptedTracks(lESDevent);

  }//if (fAnalysisType == "ESD")
  
  
  if(fAnalysisType == "AOD"){  
	  lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() ); 
	if (!lAODevent) {
		AliWarning("ERROR: lAODevent not available \n");
		return;
	}
	ncascades                      = lAODevent->GetNumberOfCascades();
        nTrackWithTPCrefitMultiplicity = -1;
        nTrackPrimaryMultiplicity      = -1;
        
        fHistCascadeMultiplicityBeforeTrigSel->Fill ( ncascades );
  }

  // For AOD or ESD ...
        nTrackMultiplicity = (InputEvent())->GetNumberOfTracks();


        //-------------------------------------------------
        fHistTrackMultiplicityForTrigEvt         ->Fill( nTrackMultiplicity             );
        fHistTPCrefitTrackMultiplicityForTrigEvt ->Fill( nTrackWithTPCrefitMultiplicity );
        fHistPrimaryTrackMultiplicityForTrigEvt  ->Fill( nTrackPrimaryMultiplicity      );
        fHistCascadeMultiplicityForTrigEvt       ->Fill( ncascades                      );




  // ---------------------------------------------------------------
  // I - Global characteristics of the events + general histos (filled for any selected events and/or for the analysed events)

                // - I.Step 1 : Characteristics of the event : prim. Vtx + magnetic field (ESD)
                //-------------

   if(fAnalysisType == "ESD"){
        const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
                // get the vtx stored in ESD found with tracks
                lPrimaryTrackingESDVtx->GetXYZ( lTrkgPrimaryVtxPos );
        
        
        const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex();	
                // get the best primary vertex available for the event
                // As done in AliCascadeVertexer, we keep the one which is the best one available.
                // between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
                // This one will be used for next calculations (DCA essentially)
                lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
                lStatusTrackingPrimVtx  = lPrimaryTrackingESDVtx->GetStatus();

        // FIXME : quality cut on the z-position of the prim vertex.
        if(fkQualityCutZprimVtxPos) {
                if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0 ) { 
                        AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !"); 
                        PostData(1, fListHistCascade); 
                        return; 
                }
        }
        
        fHistCascadeMultiplicityForTrigEvtAndZprimVtx->Fill( ncascades  );
        
        // FIXME : quality selection regarding pile-up rejection 
        if(fkRejectEventPileUp) {
                if(lESDevent->IsPileupFromSPD() ){// minContributors=3, minZdist=0.8, nSigmaZdist=3., nSigmaDiamXY=2., nSigmaDiamZ=5.  -> see http://alisoft.cern.ch/viewvc/trunk/STEER/AliESDEvent.h?root=AliRoot&r1=41914&r2=42199&pathrev=42199
                        AliWarning("Pb / Event tagged as pile-up by SPD... return !"); 
                        PostData(1, fListHistCascade); 
                        return; 
                }
        }
        
        fHistCascadeMultiplicityForTrigEvtNonPiledUpAndZprimVtx->Fill( ncascades  );
        
        // FIXME : remove TPC-only primary vertex : retain only events with tracking + SPD vertex
        if(fkQualityCutNoTPConlyPrimVtx) {
                const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
                if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtx->GetStatus() ){
                        AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
                        PostData(1, fListHistCascade); 
                        return;
                }
        }
        
        // NOTE : For older evts
        
        // As previously done in AliCascadeVertexer, we keep, between both retrieved vertices (SPD or Tracking) 
        // the one which is the best one available.
        // This one will be used for next calculations (DCA essentially)
        // At that time, the TPC-only primary vertex was not considered
        
        
        lMagneticField = lESDevent->GetMagneticField( );
        // FIXME if(TMath::Abs(lMagneticField ) < 10e-6) continue;
        
   }// end if(ESD)
        
   if(fAnalysisType == "AOD"){
        // To be developed
        const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();	
	// get the best primary vertex available for the event
	// We may keep the one which is the best one available = GetVertex(0)
	// Pb with pile-up to expect
	// This one will be used for next calculations (DCA essentially)
	lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );
        
	lStatusTrackingPrimVtx  = -1;
        lTrkgPrimaryVtxPos[0]   = -100.0;
	lTrkgPrimaryVtxPos[1]   = -100.0;
	lTrkgPrimaryVtxPos[2]   = -100.0;
        lMagneticField = 0.;
   }


                // - I.Step 2 : Filling histos that characterize the selected event : x,y,z prim. Vtx distrib. (ESD)
                //-------------

        fHistCascadeMultiplicityForSelEvt ->Fill( ncascades );
        fHistPosBestPrimaryVtxXForSelEvt  ->Fill( lBestPrimaryVtxPos[0] );
        fHistPosBestPrimaryVtxYForSelEvt  ->Fill( lBestPrimaryVtxPos[1] );
        fHistPosBestPrimaryVtxZForSelEvt  ->Fill( lBestPrimaryVtxPos[2] );
       

  
  // ---------------------------------------------------------------
  // II - Calcultaion Part dedicated to Xi vertices
  
  for (Int_t iXi = 0; iXi < ncascades; iXi++)
  {// This is the begining of the Cascade loop (ESD or AOD)
	   
    // -------------------------------------
    // II.Init - Initialisation of the local variables that will be needed for ESD/AOD

  
        // - 0th part of initialisation : around primary vertex ...
	
	Double_t lTrkgPrimaryVtxRadius3D = -500.0;
	Double_t lBestPrimaryVtxRadius3D = -500.0;

        // - 1st part of initialisation : variables needed to store AliESDCascade data members
	Double_t lEffMassXi      = 0. ;
	Double_t lChi2Xi         = -1. ;
	Double_t lDcaXiDaughters = -1. ;
	Double_t lXiCosineOfPointingAngle = -1. ;
	Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
	Double_t lXiRadius2D = -1000. ;
        Double_t lXiRadius3D = -1000. ;
        
        // - 2nd part of initialisation : Nbr of clusters within TPC for the 3 daughter cascade tracks
        Int_t    lPosTPCClusters    = -1; // For ESD only ...//FIXME : wait for availability in AOD
        Int_t    lNegTPCClusters    = -1; // For ESD only ...
        Int_t    lBachTPCClusters   = -1; // For ESD only ...
        
        Double_t lInnerWallMomCascDghters[3] = {-100., -100., -100.};
        Double_t lTPCSignalCascDghters   [3] = {-100., -100., -100.};
        
        
        // - 3rd part of initialisation : about V0 part in cascades
	Double_t lInvMassLambdaAsCascDghter = 0.;
	Double_t lV0Chi2Xi         = -1. ;
	Double_t lDcaV0DaughtersXi = -1.;
		
	Double_t lDcaBachToPrimVertexXi = -1., lDcaV0ToPrimVertexXi = -1.;
	Double_t lDcaPosToPrimVertexXi  = -1.;
	Double_t lDcaNegToPrimVertexXi  = -1.;
	Double_t lV0CosineOfPointingAngleXi = -1. ;
	Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
	Double_t lV0RadiusXi = -1000.0;
	Double_t lV0quality  = 0.;

	
	// - 4th part of initialisation : Effective masses
	Double_t lInvMassXiMinus    = 0.;
	Double_t lInvMassXiPlus     = 0.;
	Double_t lInvMassOmegaMinus = 0.;
	Double_t lInvMassOmegaPlus  = 0.;
  
	// - 5th part of initialisation : PID treatment
	Bool_t   lIsPosInXiProton      = kFALSE;
	Bool_t   lIsPosInXiPion        = kFALSE;
	Bool_t   lIsPosInOmegaProton   = kFALSE;
	Bool_t   lIsPosInOmegaPion     = kFALSE;
			
	Bool_t   lIsNegInXiProton      = kFALSE;
	Bool_t   lIsNegInXiPion        = kFALSE;
	Bool_t   lIsNegInOmegaProton   = kFALSE;
	Bool_t   lIsNegInOmegaPion     = kFALSE;
	
	Bool_t   lIsBachelorKaon       = kFALSE;
	Bool_t   lIsBachelorPion       = kFALSE; 
	
	Bool_t   lIsBachelorKaonForTPC = kFALSE; // For ESD only ...//FIXME : wait for availability in AOD
	Bool_t   lIsBachelorPionForTPC = kFALSE; // For ESD only ...
	Bool_t   lIsNegPionForTPC      = kFALSE; // For ESD only ...
	Bool_t   lIsPosPionForTPC      = kFALSE; // For ESD only ...
	Bool_t   lIsNegProtonForTPC    = kFALSE; // For ESD only ...
	Bool_t   lIsPosProtonForTPC    = kFALSE; // For ESD only ...

	// - 6th part of initialisation : extra info for QA
	Double_t lXiMomX       = 0. , lXiMomY = 0., lXiMomZ = 0.;
	Double_t lXiTransvMom  = 0. ;
	Double_t lXiTotMom     = 0. ;
		
	Double_t lBachMomX       = 0., lBachMomY  = 0., lBachMomZ   = 0.;
	Double_t lBachTransvMom  = 0.;
	Double_t lBachTotMom     = 0.;
	
	Short_t  lChargeXi = -2;
	Double_t lV0toXiCosineOfPointingAngle = 0. ;
	
	Double_t lRapXi   = -20.0, lRapOmega = -20.0,  lEta = -20.0, lTheta = 360., lPhi = 720. ;
	Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
	
  	// - 7th part of initialisation : variables for the AliCFContainer dedicated to cascade cut optmisiation
	Int_t    lSPDTrackletsMultiplicity = -1;
  
  	// - 8th part of initialisation : variables needed for Angular correlations
	TVector3 lTVect3MomXi(0.,0.,0.);
	Int_t    lArrTrackID[3] = {-1, -1, -1};

	  
  if(fAnalysisType == "ESD"){ 
  
  // -------------------------------------
  // II.ESD - Calcultaion Part dedicated to Xi vertices (ESD)
  
	AliESDcascade *xi = lESDevent->GetCascade(iXi);
	if (!xi) continue;
        
        
                // - II.Step 1 : around primary vertex
                //-------------
        lTrkgPrimaryVtxRadius3D = TMath::Sqrt(  lTrkgPrimaryVtxPos[0] * lTrkgPrimaryVtxPos[0] +
                                                lTrkgPrimaryVtxPos[1] * lTrkgPrimaryVtxPos[1] +
                                                lTrkgPrimaryVtxPos[2] * lTrkgPrimaryVtxPos[2] );

        lBestPrimaryVtxRadius3D = TMath::Sqrt(  lBestPrimaryVtxPos[0] * lBestPrimaryVtxPos[0] +
                                                lBestPrimaryVtxPos[1] * lBestPrimaryVtxPos[1] +
                                                lBestPrimaryVtxPos[2] * lBestPrimaryVtxPos[2] );


	
		// - II.Step 2 : Assigning the necessary variables for specific AliESDcascade data members (ESD)	
		//-------------
	lV0quality = 0.;
	xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis : cascade = Xi- decay

	lEffMassXi  			= xi->GetEffMassXi();
	lChi2Xi 			= xi->GetChi2Xi();
	lDcaXiDaughters			= xi->GetDcaXiDaughters();
	lXiCosineOfPointingAngle 	= xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                               lBestPrimaryVtxPos[1],
                                                                               lBestPrimaryVtxPos[2] );
		// Take care : the best available vertex should be used (like in AliCascadeVertexer)
	
	xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] ); 
	lXiRadius2D    = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
        lXiRadius3D    = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] +  lPosXi[2]*lPosXi[2]);
		
		

		// - II.Step 3 : around the tracks : Bach + V0 (ESD)
		// ~ Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
		//-------------
		
        UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
        UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
        UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
                // Care track label can be negative in MC production (linked with the track quality)
                // However = normally, not the case for track index ...
        
                // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
                if(lBachIdx == lIdxNegXi) {
                        AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
                }
                if(lBachIdx == lIdxPosXi) {
                        AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
                }
        
	AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
	AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
	AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx );
                if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
                        AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
                        continue;
                }
        
        
        lPosTPCClusters   = pTrackXi->GetTPCNcls();
        lNegTPCClusters   = nTrackXi->GetTPCNcls();
        lBachTPCClusters  = bachTrackXi->GetTPCNcls();
        
                // FIXME : rejection of a poor quality tracks
        if(fkQualityCutTPCrefit){
                // 1 - Poor quality related to TPCrefit
                ULong_t pStatus    = pTrackXi->GetStatus();
                ULong_t nStatus    = nTrackXi->GetStatus();
                ULong_t bachStatus = bachTrackXi->GetStatus();
                if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
                if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
                if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
        }
        if(fkQualityCut80TPCcls){
                // 2 - Poor quality related to TPC clusters
                if(lPosTPCClusters  < 80) { AliWarning("Pb / V0 Pos. track has less than 80 TPC clusters ... continue!"); continue; }
                if(lNegTPCClusters  < 80) { AliWarning("Pb / V0 Neg. track has less than 80 TPC clusters ... continue!"); continue; }
                if(lBachTPCClusters < 80) { AliWarning("Pb / Bach.   track has less than 80 TPC clusters ... continue!"); continue; }
        }
        
        const AliExternalTrackParam *pExtTrack    = pTrackXi    ->GetInnerParam();
        const AliExternalTrackParam *nExtTrack    = nTrackXi    ->GetInnerParam();
        const AliExternalTrackParam *bachExtTrack = bachTrackXi ->GetInnerParam();
        
        if (pExtTrack) {
                lInnerWallMomCascDghters[0] = pExtTrack ->GetP() * pExtTrack ->Charge();
                lTPCSignalCascDghters   [0] = pTrackXi  ->GetTPCsignal();
        }
        if (nExtTrack) {
                lInnerWallMomCascDghters[1] = nExtTrack ->GetP() * nExtTrack ->Charge();
                lTPCSignalCascDghters   [1] = nTrackXi  ->GetTPCsignal();
        }
	if (bachExtTrack) {
                lInnerWallMomCascDghters[2] = bachExtTrack ->GetP() * bachExtTrack ->Charge();
                lTPCSignalCascDghters   [2] = bachTrackXi  ->GetTPCsignal();
        }


	lInvMassLambdaAsCascDghter	= xi->GetEffMass();
		// This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
	lDcaV0DaughtersXi 		= xi->GetDcaV0Daughters(); 
	lV0Chi2Xi 			= xi->GetChi2V0();
	
	lV0CosineOfPointingAngleXi 	= xi->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0],
									  lBestPrimaryVtxPos[1],
									  lBestPrimaryVtxPos[2] );

	lDcaV0ToPrimVertexXi 		= xi->GetD( lBestPrimaryVtxPos[0], 
						    lBestPrimaryVtxPos[1], 
						    lBestPrimaryVtxPos[2] );
		
	lDcaBachToPrimVertexXi = TMath::Abs( bachTrackXi->GetD(	lBestPrimaryVtxPos[0], 
						     		lBestPrimaryVtxPos[1], 
						     		lMagneticField  ) ); 
					// Note : AliExternalTrackParam::GetD returns an algebraic value ...
		
		xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] ); 
	lV0RadiusXi		= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
	
	lDcaPosToPrimVertexXi 	= TMath::Abs( pTrackXi	->GetD(	lBestPrimaryVtxPos[0], 
						 		lBestPrimaryVtxPos[1], 
						 		lMagneticField  )     ); 
	
	lDcaNegToPrimVertexXi 	= TMath::Abs( nTrackXi	->GetD(	lBestPrimaryVtxPos[0], 
					    			lBestPrimaryVtxPos[1], 
					    			lMagneticField  )     ); 
		
		// - II.Step 3' : extra-selection for cascade candidates
		// Towards optimisation of AA selection
        // FIXME
	if(fkExtraSelections){
		// if(lChi2Xi > 2000) continue;
		// if(lV0Chi2Xi > 2000) continue;
		
		if(lDcaXiDaughters > 0.05) continue; // > 0.1 by default
		//if(lXiCosineOfPointingAngle < 0.999 ) continue;
		if(lXiRadius2D < 1.0) continue;
		if(lXiRadius2D > 100) continue;
		if(TMath::Abs(lInvMassLambdaAsCascDghter-1.11568) > 0.008) continue;
		if(lDcaV0DaughtersXi > 0.3) continue;
		
		if(lV0CosineOfPointingAngleXi > 0.9999) continue;
		//if(lDcaV0ToPrimVertexXi < 0.09) continue;
		if(lDcaBachToPrimVertexXi < 0.04) continue;
		
		if(lV0RadiusXi < 1.0) continue;
		if(lV0RadiusXi > 100) continue;
		//if(lDcaPosToPrimVertexXi < 0.6) continue;
		//if(lDcaNegToPrimVertexXi < 0.6) continue;
	}
	
	
	
		// - II.Step 4 : around effective masses (ESD)
		// ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
		//-------------

	
	if( bachTrackXi->Charge() < 0 )	{
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , 3312); 	
			// Calculate the effective mass of the Xi- candidate. 
			// pdg code 3312 = Xi-
		lInvMassXiMinus = xi->GetEffMassXi();
		
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , 3334); 	
			// Calculate the effective mass of the Xi- candidate. 
			// pdg code 3334 = Omega-
		lInvMassOmegaMinus = xi->GetEffMassXi();
					
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , 3312); 	// Back to default hyp.
	}// end if negative bachelor
	
	
	if( bachTrackXi->Charge() >  0 ){
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , -3312); 	
			// Calculate the effective mass of the Xi+ candidate. 
			// pdg code -3312 = Xi+
		lInvMassXiPlus = xi->GetEffMassXi();
		
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , -3334); 	
			// Calculate the effective mass of the Xi+ candidate. 
			// pdg code -3334  = Omega+
		lInvMassOmegaPlus = xi->GetEffMassXi();
		
		lV0quality = 0.;
		xi->ChangeMassHypothesis(lV0quality , -3312); 	// Back to "default" hyp.
	}// end if positive bachelor
	
	
	
		// - II.Step 5 : PID on the daughter tracks
		//-------------
	
	// A - Combined PID
	// Reasonable guess for the priors for the cascade track sample (e-, mu, pi, K, p)// Bo: neutral part added now !
	Double_t lPriorsGuessXi[10]    = {0, 0, 2, 0, 1, 0, 0, 0, 0, 0};
	Double_t lPriorsGuessOmega[10] = {0, 0, 1, 1, 1, 0, 0, 0, 0, 0};
	
	// Combined VO-positive-daughter PID
	AliPID pPidXi;		pPidXi.SetPriors(    lPriorsGuessXi,    1); // Bo: now needed to specify charged
	AliPID pPidOmega;	pPidOmega.SetPriors( lPriorsGuessOmega, 1); // Bo: now needed to specify charged
		
	if( pTrackXi->IsOn(AliESDtrack::kESDpid) ){  // Combined PID exists
		Double_t r[10] = {0.}; pTrackXi->GetESDpid(r);
		pPidXi.SetProbabilities(r);
		pPidOmega.SetProbabilities(r);
		
		// Check if the V0 positive track is a proton (case for Xi-)
		Double_t pproton = pPidXi.GetProbability(AliPID::kProton);
		if (pproton > pPidXi.GetProbability(AliPID::kElectron) &&
		    pproton > pPidXi.GetProbability(AliPID::kMuon)     &&
		    pproton > pPidXi.GetProbability(AliPID::kPion)     &&
		    pproton > pPidXi.GetProbability(AliPID::kKaon)     )     lIsPosInXiProton = kTRUE;
		
		// Check if the V0 positive track is a pi+ (case for Xi+)
		Double_t ppion = pPidXi.GetProbability(AliPID::kPion);
		if (ppion > pPidXi.GetProbability(AliPID::kElectron) &&
		    ppion > pPidXi.GetProbability(AliPID::kMuon)     &&
		    ppion > pPidXi.GetProbability(AliPID::kKaon)     &&
		    ppion > pPidXi.GetProbability(AliPID::kProton)   )     lIsPosInXiPion = kTRUE;
		
		
		// Check if the V0 positive track is a proton (case for Omega-)
		pproton = 0.;
		    pproton = pPidOmega.GetProbability(AliPID::kProton);
		if (pproton > pPidOmega.GetProbability(AliPID::kElectron) &&
		    pproton > pPidOmega.GetProbability(AliPID::kMuon)     &&
		    pproton > pPidOmega.GetProbability(AliPID::kPion)     &&
		    pproton > pPidOmega.GetProbability(AliPID::kKaon)     )  lIsPosInOmegaProton = kTRUE;
		
		// Check if the V0 positive track is a pi+ (case for Omega+)
		ppion = 0.;
		    ppion = pPidOmega.GetProbability(AliPID::kPion);
		if (ppion > pPidOmega.GetProbability(AliPID::kElectron) &&
		    ppion > pPidOmega.GetProbability(AliPID::kMuon)     &&
		    ppion > pPidOmega.GetProbability(AliPID::kKaon)     &&
		    ppion > pPidOmega.GetProbability(AliPID::kProton)   )    lIsPosInOmegaPion = kTRUE;
		
	}// end if V0 positive track with existing combined PID	
	
	
	// Combined VO-negative-daughter PID
	AliPID nPidXi;		nPidXi.SetPriors(    lPriorsGuessXi,    1); // Bo: now needed to specify charged
	AliPID nPidOmega;	nPidOmega.SetPriors( lPriorsGuessOmega, 1); // Bo: now needed to specify charged
		
	if( nTrackXi->IsOn(AliESDtrack::kESDpid) ){  // Combined PID exists
		Double_t r[10] = {0.}; nTrackXi->GetESDpid(r);
		nPidXi.SetProbabilities(r);
		nPidOmega.SetProbabilities(r);
		
		// Check if the V0 negative track is a pi- (case for Xi-)
		Double_t ppion = nPidXi.GetProbability(AliPID::kPion);
		if (ppion > nPidXi.GetProbability(AliPID::kElectron) &&
		    ppion > nPidXi.GetProbability(AliPID::kMuon)     &&
		    ppion > nPidXi.GetProbability(AliPID::kKaon)     &&
		    ppion > nPidXi.GetProbability(AliPID::kProton)   )     lIsNegInXiPion = kTRUE;

		// Check if the V0 negative track is an anti-proton (case for Xi+)
		Double_t pproton = nPidXi.GetProbability(AliPID::kProton);
		if (pproton > nPidXi.GetProbability(AliPID::kElectron) &&
		    pproton > nPidXi.GetProbability(AliPID::kMuon)     &&
		    pproton > nPidXi.GetProbability(AliPID::kPion)     &&
		    pproton > nPidXi.GetProbability(AliPID::kKaon)     )     lIsNegInXiProton = kTRUE;
		
		// Check if the V0 negative track is a pi- (case for Omega-)
		ppion = 0.;
		    ppion = nPidOmega.GetProbability(AliPID::kPion);
		if (ppion > nPidOmega.GetProbability(AliPID::kElectron) &&
		    ppion > nPidOmega.GetProbability(AliPID::kMuon)     &&
		    ppion > nPidOmega.GetProbability(AliPID::kKaon)     &&
		    ppion > nPidOmega.GetProbability(AliPID::kProton)   )    lIsNegInOmegaPion = kTRUE;
		
		// Check if the V0 negative track is an anti-proton (case for Omega+)
		pproton = 0.;
		         pproton = nPidOmega.GetProbability(AliPID::kProton);
		if (pproton > nPidOmega.GetProbability(AliPID::kElectron) &&
		    pproton > nPidOmega.GetProbability(AliPID::kMuon)     &&
		    pproton > nPidOmega.GetProbability(AliPID::kPion)     &&
		    pproton > nPidOmega.GetProbability(AliPID::kKaon)     )  lIsNegInOmegaProton = kTRUE;
		
	}// end if V0 negative track with existing combined PID	
	
		
	// Combined bachelor PID
	AliPID bachPidXi;	bachPidXi.SetPriors(    lPriorsGuessXi,    1); // Bo: now needed to specify charged
	AliPID bachPidOmega;	bachPidOmega.SetPriors( lPriorsGuessOmega, 1); // Bo: now needed to specify charged
	
	if( bachTrackXi->IsOn(AliESDtrack::kESDpid) ){  // Combined PID exists
		Double_t r[10] = {0.}; bachTrackXi->GetESDpid(r);
		bachPidXi.SetProbabilities(r);
		bachPidOmega.SetProbabilities(r);
		// Check if the bachelor track is a pion
		Double_t ppion = bachPidXi.GetProbability(AliPID::kPion);
		if (ppion > bachPidXi.GetProbability(AliPID::kElectron) &&
		    ppion > bachPidXi.GetProbability(AliPID::kMuon)     &&
		    ppion > bachPidXi.GetProbability(AliPID::kKaon)     &&
		    ppion > bachPidXi.GetProbability(AliPID::kProton)   )     lIsBachelorPion = kTRUE;
		// Check if the bachelor track is a kaon
		Double_t pkaon = bachPidOmega.GetProbability(AliPID::kKaon);
		if (pkaon > bachPidOmega.GetProbability(AliPID::kElectron) &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kMuon)     &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kPion)     &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kProton)   )  lIsBachelorKaon = kTRUE;	
	}// end if bachelor track with existing combined PID
	
	
	// B - TPC PID : 3-sigma bands on Bethe-Bloch curve
	
        // Bachelor
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 4) lIsBachelorKaonForTPC = kTRUE;
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 4) lIsBachelorPionForTPC = kTRUE;
        
        // Negative V0 daughter
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 4) lIsNegPionForTPC   = kTRUE;
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 4) lIsNegProtonForTPC = kTRUE;
        
        // Positive V0 daughter
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 4) lIsPosPionForTPC   = kTRUE;
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 4) lIsPosProtonForTPC = kTRUE;
        
        /*
        const AliExternalTrackParam *pInnerWallTrackXi    = pTrackXi    ->GetInnerParam();
        const AliExternalTrackParam *nInnerWallTrackXi    = nTrackXi    ->GetInnerParam();
        const AliExternalTrackParam *bachInnerWallTrackXi = bachTrackXi ->GetInnerParam();
        if(pInnerWallTrackXi && nInnerWallTrackXi && bachInnerWallTrackXi ){
                
                Double_t pMomInnerWall    = pInnerWallTrackXi   ->GetP();
                Double_t nMomInnerWall    = nInnerWallTrackXi   ->GetP();
                Double_t bachMomInnerWall = bachInnerWallTrackXi->GetP();
                
                // Bachelor
                if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 3)                              lIsBachelorPionForTPC = kTRUE;
                if (bachMomInnerWall < 0.350  && TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 5) lIsBachelorKaonForTPC = kTRUE;
                if (bachMomInnerWall > 0.350  && TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 3) lIsBachelorKaonForTPC = kTRUE;
                
                // Negative V0 daughter
                if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 3  )                           lIsNegPionForTPC   = kTRUE;
                if (nMomInnerWall < 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton ) ) < 5 )   lIsNegProtonForTPC = kTRUE;
                if (nMomInnerWall > 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton ) ) < 3 )   lIsNegProtonForTPC = kTRUE;
                
                // Positive V0 daughter
                if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 3 )                            lIsPosPionForTPC   = kTRUE;
                if (pMomInnerWall < 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 5)     lIsPosProtonForTPC = kTRUE;
                if (pMomInnerWall > 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 3)     lIsPosProtonForTPC = kTRUE;
        }
        */
        
        
        
		
		// - II.Step 6 : extra info for QA (ESD)
		// miscellaneous pieces of info that may help regarding data quality assessment.
		//-------------

	xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
		lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
		lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
		
	xi->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
		lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
		lBachTotMom  	= TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );

	lChargeXi = xi->Charge();

	lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
	
	lRapXi    = xi->RapXi();
	lRapOmega = xi->RapOmega();
	lEta      = xi->Eta();
	lTheta    = xi->Theta() *180.0/TMath::Pi();
	lPhi      = xi->Phi()   *180.0/TMath::Pi();
	lAlphaXi  = xi->AlphaXi();
	lPtArmXi  = xi->PtArmXi();
	
	
	//FIXME : Extra-cut = Anti-splitting cut for lambda daughters
	Bool_t kAntiSplittingLambda = kFALSE;
	
	if(kAntiSplittingLambda){
		Double_t lNMomX = 0., lNMomY = 0., lNMomZ = 0.;
		Double_t lPMomX = 0., lPMomY = 0., lPMomZ = 0.;
		
		xi->GetPPxPyPz(lPMomX, lPMomY, lPMomZ); 
		xi->GetNPxPyPz(lNMomX, lNMomY, lNMomZ); 
		
		if( xi->Charge() < 0){// Xi- or Omega-
			if(TMath::Abs(lBachTransvMom - TMath::Sqrt( lNMomX*lNMomX + lNMomY*lNMomY )  ) < 0.075) continue;
		}
		else {                //Xi+ or Omega+
			if(TMath::Abs(lBachTransvMom - TMath::Sqrt( lPMomX*lPMomX + lPMomY*lPMomY ) ) < 0.075) continue;
		}
	}
	
        //FIXME : Just to know which file is currently open : locate the file containing Xi
        // cout << "Name of the file containing Xi candidate(s) :" 
        //        << CurrentFileName() 
        //        << " / entry: "     << Entry()
        //        << " / in file: "   << lESDevent->GetEventNumberInFile()   // <- Cvetan / From Mihaela: AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetReadEntry();
        //        << " : mass(Xi) = " << xi->GetEffMassXi() 
        //        << " / charge = "   << lChargeXi
        //        << " / pt(Casc) = " << lXiTransvMom
        //        << " / Decay 2d R(Xi) = " << lXiRadius2D 
        //        << " / Track Index(Pos)  = " << lIdxPosXi << "/ Nb(TPC clusters) = " << lPosTPCClusters 
        //        << " / Track Index(Neg)  = " << lIdxNegXi << "/ Nb(TPC clusters) = " << lNegTPCClusters 
        //        << " / Track Index(Bach) = " << lBachIdx  << "/ Nb(TPC clusters) = " << lBachTPCClusters 
        //        << endl;

	
		// II.Step 7 - Complementary info for monitoring the cascade cut variables
	
	const AliMultiplicity *lAliMult = lESDevent->GetMultiplicity();
	lSPDTrackletsMultiplicity = lAliMult->GetNumberOfTracklets();
	
        	// II.Step 8 - Azimuthal correlation study
		//-------------
	
	lTVect3MomXi.SetXYZ( lXiMomX, lXiMomY, lXiMomZ );
	lArrTrackID[0] = pTrackXi   ->GetID();
	lArrTrackID[1] = nTrackXi   ->GetID();
	lArrTrackID[2] = bachTrackXi->GetID();
	
	
  }// end of ESD treatment
  
 
  if(fAnalysisType == "AOD"){
	
	// -------------------------------------
	// II.AOD - Calcultaion Part dedicated to Xi vertices (ESD)
	
	const AliAODcascade *xi = lAODevent->GetCascade(iXi);
	if (!xi) continue;
		
	// Just to know which file is currently open : locate the file containing Xi
	// cout << "Name of the file containing Xi candidate(s) :" <<  fesdH->GetTree()->GetCurrentFile()->GetName() << endl;
	
	
		// - II.Step 1 : Characteristics of the event : prim. Vtx + magnetic field (AOD)
		//-------------
	

	lTrkgPrimaryVtxRadius3D = -500. ;
	// FIXME : We don't have the different prim. vertex at the AOD level -> nothing to do.

	lBestPrimaryVtxRadius3D = TMath::Sqrt(  lBestPrimaryVtxPos[0] * lBestPrimaryVtxPos[0] +
						lBestPrimaryVtxPos[1] * lBestPrimaryVtxPos[1] +
						lBestPrimaryVtxPos[2] * lBestPrimaryVtxPos[2] );
		
	
		// - II.Step 2 : Assigning the necessary variables for specific AliAODcascade data members (AOD)	
		//-------------
	
	lEffMassXi  			= xi->MassXi(); // default working hypothesis : cascade = Xi- decay
	lChi2Xi 			= xi->Chi2Xi();
	lDcaXiDaughters			= xi->DcaXiDaughters();
	lXiCosineOfPointingAngle 	= xi->CosPointingAngleXi( lBestPrimaryVtxPos[0], 
								  lBestPrimaryVtxPos[1], 
								  lBestPrimaryVtxPos[2] );
					// Take care : 
					// the best available vertex should be used (like in AliCascadeVertexer)

		lPosXi[0] = xi->DecayVertexXiX();
		lPosXi[1] = xi->DecayVertexXiY();
		lPosXi[2] = xi->DecayVertexXiZ();
	lXiRadius2D = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
        lXiRadius3D = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] +  lPosXi[2]*lPosXi[2] );
		

		// - II.Step 3 : around the tracks : Bach + V0 (AOD)
		// ~ Necessary variables for AODcascade data members coming from the AODv0 part (inheritance)
		//-------------
		
	lChargeXi 			= xi->ChargeXi();
	
	if( lChargeXi < 0)  	
	  lInvMassLambdaAsCascDghter	= xi->MassLambda();
	else 			
	  lInvMassLambdaAsCascDghter	= xi->MassAntiLambda();

	lDcaV0DaughtersXi 		= xi->DcaV0Daughters(); 
	lV0Chi2Xi 			= xi->Chi2V0();

	lV0CosineOfPointingAngleXi 	= xi->CosPointingAngle( lBestPrimaryVtxPos );
	lDcaV0ToPrimVertexXi 		= xi->DcaV0ToPrimVertex();
	
	lDcaBachToPrimVertexXi 		= xi->DcaBachToPrimVertex(); 
	
	
		lPosV0Xi[0] = xi->DecayVertexV0X();
		lPosV0Xi[1] = xi->DecayVertexV0Y();
		lPosV0Xi[2] = xi->DecayVertexV0Z(); 
	lV0RadiusXi	= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );

	lDcaPosToPrimVertexXi		= xi->DcaPosToPrimVertex(); 
	lDcaNegToPrimVertexXi		= xi->DcaNegToPrimVertex(); 


		// - II.Step 4 : around effective masses (AOD)
		// ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
		//-------------

	if( lChargeXi < 0 )		lInvMassXiMinus 	= xi->MassXi();
	if( lChargeXi > 0 )		lInvMassXiPlus 		= xi->MassXi();
	if( lChargeXi < 0 )		lInvMassOmegaMinus 	= xi->MassOmega();
	if( lChargeXi > 0 )		lInvMassOmegaPlus 	= xi->MassOmega();

	
		// - II.Step 5 : PID on the daughters (To be developed ...)
		//-------------
	
	// Combined PID
	
	/* 
	// Reasonable guess for the priors for the cascade track sample
	Double_t lPriorsGuessXi[5]    = {0.0, 0.0, 2, 0, 1};
	Double_t lPriorsGuessOmega[5] = {0.0, 0.0, 1, 1, 1};
	AliPID bachPidXi;	bachPidXi.SetPriors(    lPriorsGuessXi    );
	AliPID bachPidOmega;	bachPidOmega.SetPriors( lPriorsGuessOmega );
	
	const AliAODTrack *bachTrackXi = lAODevent->GetTrack( xi->GetBachID() ); // FIXME : GetBachID not implemented ?
	
	if( bachTrackXi->IsOn(AliESDtrack::kESDpid) ){  // Combined PID exists, the AOD flags = a copy of the ESD ones
		Double_t r[10]; bachTrackXi->GetPID(r);
		bachPidXi.SetProbabilities(r);
		bachPidOmega.SetProbabilities(r);
		// Check if the bachelor track is a pion
		Double_t ppion = bachPidXi.GetProbability(AliPID::kPion);
		if (ppion > bachPidXi.GetProbability(AliPID::kElectron) &&
		    ppion > bachPidXi.GetProbability(AliPID::kMuon)     &&
		    ppion > bachPidXi.GetProbability(AliPID::kKaon)     &&
		    ppion > bachPidXi.GetProbability(AliPID::kProton)   )     lIsBachelorPion = kTRUE;
		// Check if the bachelor track is a kaon
		Double_t pkaon = bachPidOmega.GetProbability(AliPID::kKaon);
		if (pkaon > bachPidOmega.GetProbability(AliPID::kElectron) &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kMuon)     &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kPion)     &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kProton)   )  lIsBachelorKaon = kTRUE;
		
	}// end if bachelor track with existing combined PID
	*/
	
	// TPC PID
	
		// - II.Step 6 : extra info for QA (AOD)
		// miscellaneous pieces onf info that may help regarding data quality assessment.
		//-------------

		lXiMomX = xi->MomXiX();
		lXiMomY = xi->MomXiY();
		lXiMomZ = xi->MomXiZ();
	lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
	lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
	
	 	lBachMomX = xi->MomBachX();
	 	lBachMomY = xi->MomBachY();
	 	lBachMomZ = xi->MomBachZ();		
	lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
	lBachTotMom  	= TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );

	
	lV0toXiCosineOfPointingAngle = xi->CosPointingAngle( xi->GetDecayVertexXi() );
	
	lRapXi    = xi->RapXi();
	lRapOmega = xi->RapOmega();
	lEta      = xi->Eta();				// Will not work ! need a method Pz(), Py() Px() 
	lTheta    = xi->Theta() *180.0/TMath::Pi();     // in AODcascade.
	lPhi      = xi->Phi()   *180.0/TMath::Pi();     // Here, we will get eta, theta, phi for the V0 ...
	lAlphaXi  = xi->AlphaXi();
	lPtArmXi  = xi->PtArmXi();

		// II.Step 7 - Complementary info for monitoring the cascade cut variables
	//FIXME : missing for AOD : Tacklet Multiplicity + TPCCluster
	
		// II.Step 8 - Azimuthal correlation study
		//-------------
	
	lTVect3MomXi.SetXYZ( lXiMomX, lXiMomY, lXiMomZ );
	
	AliAODTrack *pTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
	AliAODTrack *nTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
	AliAODTrack *bachTrackXi = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );	
		if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
			AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
			continue;
		}
		
	lArrTrackID[0] = pTrackXi   ->GetID();
	lArrTrackID[1] = nTrackXi   ->GetID();
	lArrTrackID[2] = bachTrackXi->GetID();
	
  }// end of AOD treatment


    // -------------------------------------
    // II.Fill - Filling the TH1,2,3Fs, HnSparses, CFContainers, FOR events with CASCADES !
        
        // if( lIsBachelorKaonForTPC )
        //                 // FIXME : Just to know which file is currently open : locate the file containing Xi
        // cout << "Name of the file containing Omega candidate(s) :" 
        //         << CurrentFileName() 
        //         << " / entry: "     << Entry()
        //         << " / in file: "   << lESDevent->GetEventNumberInFile()   // <- Cvetan / From Mihaela: AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetReadEntry();
        //         << " : mass(Omega+) = " << lInvMassOmegaPlus 
        //         << " : mass(Omega-) = " << lInvMassOmegaMinus
        //         << " / charge = "   << lChargeXi
        //         << " / pt(Casc) = " << lXiTransvMom
        //         << " / Decay 2d R(Xi) = " << lXiRadius2D 
        //         << endl;
  

	// - II.Fill.Step 1	 : primary vertex
  
        fHistTPCrefitTrackMultiplicityForCascadeEvt->Fill( nTrackWithTPCrefitMultiplicity );
        fHistPrimaryTrackMultiplicityForCascadeEvt ->Fill( nTrackPrimaryMultiplicity );
        
        fHistPosV0TPCClusters           ->Fill( lPosTPCClusters  );
        fHistNegV0TPCClusters           ->Fill( lNegTPCClusters  );
        fHistBachTPCClusters            ->Fill( lBachTPCClusters );
        
        f2dHistTPCdEdxOfCascDghters     ->Fill( lInnerWallMomCascDghters[0] , lTPCSignalCascDghters[0]  );
        f2dHistTPCdEdxOfCascDghters     ->Fill( lInnerWallMomCascDghters[1] , lTPCSignalCascDghters[1]  );
        f2dHistTPCdEdxOfCascDghters     ->Fill( lInnerWallMomCascDghters[2] , lTPCSignalCascDghters[2]  );
        
	fHistVtxStatus			->Fill( lStatusTrackingPrimVtx   );  // 1 if tracking vtx = ok

	if( lStatusTrackingPrimVtx ){
		fHistPosTrkgPrimaryVtxXForCascadeEvt  ->Fill( lTrkgPrimaryVtxPos[0]    );
		fHistPosTrkgPrimaryVtxYForCascadeEvt  ->Fill( lTrkgPrimaryVtxPos[1]    );
		fHistPosTrkgPrimaryVtxZForCascadeEvt  ->Fill( lTrkgPrimaryVtxPos[2]    );
		fHistTrkgPrimaryVtxRadius             ->Fill( lTrkgPrimaryVtxRadius3D );
	}

	fHistPosBestPrimaryVtxXForCascadeEvt   ->Fill( lBestPrimaryVtxPos[0]    );
	fHistPosBestPrimaryVtxYForCascadeEvt   ->Fill( lBestPrimaryVtxPos[1]    );
	fHistPosBestPrimaryVtxZForCascadeEvt   ->Fill( lBestPrimaryVtxPos[2]    );
	fHistBestPrimaryVtxRadius              ->Fill( lBestPrimaryVtxRadius3D  );
	
	f2dHistTrkgPrimVtxVsBestPrimVtx->Fill( lTrkgPrimaryVtxRadius3D, lBestPrimaryVtxRadius3D );

        // **************************** With PID on ? ... for the signal region ? ************FIXME**************************************
        if( ( (lChargeXi<0) && lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) ||
            ( (lChargeXi>0) && lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )  )
                // NOTE : 
                // with this PID condition, it could happen that a cascade candidate satisfies the wrong requirement,
                // e.g. one looks at a Xi- candidate for which lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC = kFALSE
                //      Expectation: it should be excluded.
                //      but lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC = kTRUE
                //      then this bad Xi-candidate will contribute anyway (OR condition).
                // Hence : the extra condition on the sign of the Cascade
        {
                //         if( TMath::Abs( lInvMassXiMinus-1.3217 ) < 0.010 || TMath::Abs( lInvMassXiPlus-1.3217 ) < 0.010){}  }
                        
                // II.Fill.Step 2
                fHistEffMassXi			->Fill( lEffMassXi               );
                fHistChi2Xi			->Fill( lChi2Xi                  );	// Flag CascadeVtxer: Cut Variable a
                fHistDcaXiDaughters		->Fill( lDcaXiDaughters          );	// Flag CascadeVtxer: Cut Variable e 
                fHistDcaBachToPrimVertex	->Fill( lDcaBachToPrimVertexXi   );	// Flag CascadeVtxer: Cut Variable d
                fHistXiCosineOfPointingAngle	->Fill( lXiCosineOfPointingAngle );	// Flag CascadeVtxer: Cut Variable f
                fHistXiRadius			->Fill( lXiRadius2D                );	// Flag CascadeVtxer: Cut Variable g+h
                
                
                // II.Fill.Step 3
                fHistMassLambdaAsCascDghter	->Fill( lInvMassLambdaAsCascDghter );	// Flag CascadeVtxer: Cut Variable c
                fHistV0Chi2Xi			->Fill( lV0Chi2Xi                  );	
                fHistDcaV0DaughtersXi		->Fill( lDcaV0DaughtersXi          );
                fHistV0CosineOfPointingAngleXi	->Fill( lV0CosineOfPointingAngleXi ); 
                fHistV0RadiusXi			->Fill( lV0RadiusXi                );
                
                fHistDcaV0ToPrimVertexXi	->Fill( lDcaV0ToPrimVertexXi       );	// Flag CascadeVtxer: Cut Variable b
                fHistDcaPosToPrimVertexXi	->Fill( lDcaPosToPrimVertexXi      );
                fHistDcaNegToPrimVertexXi	->Fill( lDcaNegToPrimVertexXi      );
                
        
                // II.Fill.Step 4 : extra QA info
                
                fHistChargeXi                   ->Fill( lChargeXi      );
                fHistV0toXiCosineOfPointingAngle->Fill( lV0toXiCosineOfPointingAngle );
        
                if( TMath::Abs( lInvMassXiMinus-1.3217 ) < 0.010 || TMath::Abs( lInvMassXiPlus-1.3217 ) < 0.010){// One InvMass should be different from 0
                        fHistXiTransvMom        ->Fill( lXiTransvMom   );
                        fHistXiTotMom           ->Fill( lXiTotMom      );
        
                        fHistBachTransvMomXi    ->Fill( lBachTransvMom );
                        fHistBachTotMomXi       ->Fill( lBachTotMom    );
        
                        fHistRapXi              ->Fill( lRapXi         );
                        fHistEtaXi              ->Fill( lEta           );
                        fHistThetaXi            ->Fill( lTheta         );
                        fHistPhiXi              ->Fill( lPhi           );

                }
                
                if( (lChargeXi < 0) && (TMath::Abs( lInvMassXiMinus-1.3217 ) < 0.010) )         fHistcTauXiMinus     ->Fill( lXiRadius3D * 1.3271/lXiTotMom );
                if( (lChargeXi > 0) && (TMath::Abs( lInvMassXiPlus -1.3217 ) < 0.010) )         fHistcTauXiPlus      ->Fill( lXiRadius3D * 1.3271/lXiTotMom );

                if( TMath::Abs( lInvMassOmegaMinus-1.672 ) < 0.010 || TMath::Abs( lInvMassOmegaPlus-1.672 ) < 0.010 ){// One InvMass should be different from 0
                        fHistRapOmega           ->Fill( lRapOmega            ); 

                }
                
                if( (lChargeXi < 0) && (TMath::Abs( lInvMassOmegaMinus-1.672 ) < 0.010) )        fHistcTauOmegaMinus ->Fill( lXiRadius3D * 1.67245/lXiTotMom );
                if( (lChargeXi > 0) && (TMath::Abs( lInvMassOmegaPlus- 1.672 ) < 0.010) )        fHistcTauOmegaPlus  ->Fill( lXiRadius3D * 1.67245/lXiTotMom );
        
                f2dHistArmenteros	        ->Fill( lChargeXi*lAlphaXi, lPtArmXi   );
        }// end with PID ...
	
	// II.Fill.Step 5 : inv mass plots 1D
	if( lChargeXi < 0 ){
					fHistMassXiMinus	       ->Fill( lInvMassXiMinus    );
					fHistMassOmegaMinus	       ->Fill( lInvMassOmegaMinus );
		if(lIsBachelorPion)	fHistMassWithCombPIDXiMinus    ->Fill( lInvMassXiMinus    );
		if(lIsBachelorKaon)	fHistMassWithCombPIDOmegaMinus ->Fill( lInvMassOmegaMinus );
	}
	
	if( lChargeXi > 0 ){
					fHistMassXiPlus		       ->Fill( lInvMassXiPlus     );
					fHistMassOmegaPlus	       ->Fill( lInvMassOmegaPlus  );
		if(lIsBachelorPion)	fHistMassWithCombPIDXiPlus     ->Fill( lInvMassXiPlus     );
		if(lIsBachelorKaon)	fHistMassWithCombPIDOmegaPlus  ->Fill( lInvMassOmegaPlus  );
	}
	
	
	// II.Fill.Step 6 : inv mass plots 2D, 3D
	if( lChargeXi < 0 ) {
		f2dHistEffMassLambdaVsEffMassXiMinus->Fill( lInvMassLambdaAsCascDghter, lInvMassXiMinus ); 
		f2dHistEffMassXiVsEffMassOmegaMinus ->Fill( lInvMassXiMinus, lInvMassOmegaMinus );
		f2dHistXiRadiusVsEffMassXiMinus     ->Fill( lXiRadius2D, lInvMassXiMinus );
		f2dHistXiRadiusVsEffMassOmegaMinus  ->Fill( lXiRadius2D, lInvMassOmegaMinus );
		f3dHistXiPtVsEffMassVsYXiMinus      ->Fill( lXiTransvMom, lInvMassXiMinus,    lRapXi    );
		f3dHistXiPtVsEffMassVsYOmegaMinus   ->Fill( lXiTransvMom, lInvMassOmegaMinus, lRapOmega );
	}
	else{
		f2dHistEffMassLambdaVsEffMassXiPlus ->Fill( lInvMassLambdaAsCascDghter, lInvMassXiPlus );
		f2dHistEffMassXiVsEffMassOmegaPlus  ->Fill( lInvMassXiPlus, lInvMassOmegaPlus );
		f2dHistXiRadiusVsEffMassXiPlus      ->Fill( lXiRadius2D, lInvMassXiPlus);
		f2dHistXiRadiusVsEffMassOmegaPlus   ->Fill( lXiRadius2D, lInvMassOmegaPlus );
		f3dHistXiPtVsEffMassVsYXiPlus       ->Fill( lXiTransvMom, lInvMassXiPlus,    lRapXi    );
		f3dHistXiPtVsEffMassVsYOmegaPlus    ->Fill( lXiTransvMom, lInvMassOmegaPlus, lRapOmega );
	}
	
	// - Filling the AliCFContainers related to PID
	
	Double_t lContainerPIDVars[4] = {0.0};
	
	// Xi Minus		
	if( lChargeXi < 0 ) {
		lContainerPIDVars[0] = lXiTransvMom       ;
		lContainerPIDVars[1] = lInvMassXiMinus    ;
		lContainerPIDVars[2] = lRapXi             ;
		lContainerPIDVars[3] = nTrackPrimaryMultiplicity;       // FIXME : nTrackPrimaryMultiplicity not set for AOD ... = -1
			
		// No PID
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorPionForTPC  )
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
		if( lIsBachelorPionForTPC && 
		    lIsPosProtonForTPC     )
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorPionForTPC && 
		    lIsPosProtonForTPC    && 
		    lIsNegPionForTPC       )
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorPion        )
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorPion       && 
		    lIsPosInXiProton    )
			fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if(lIsBachelorPion     && 
		   lIsPosInXiProton && 
		   lIsNegInXiPion    )
		 	fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Xi Plus		
	if( lChargeXi > 0 ) {
		lContainerPIDVars[0] = lXiTransvMom       ;
		lContainerPIDVars[1] = lInvMassXiPlus     ;
		lContainerPIDVars[2] = lRapXi             ;
		lContainerPIDVars[3] = nTrackPrimaryMultiplicity;       // FIXME : nTrackPrimaryMultiplicity not set for AOD ... = -1
			
		// No PID
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorPionForTPC  )
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
		if( lIsBachelorPionForTPC && 
		    lIsNegProtonForTPC     )
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorPionForTPC && 
		    lIsNegProtonForTPC    && 
		    lIsPosPionForTPC       )
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorPion        )
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorPion       && 
		    lIsNegInXiProton    )
			fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if(lIsBachelorPion     && 
		   lIsNegInXiProton && 
		   lIsPosInXiPion    )
		 	fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Omega Minus		
	if( lChargeXi < 0 ) {
		lContainerPIDVars[0] = lXiTransvMom       ;
		lContainerPIDVars[1] = lInvMassOmegaMinus ;
		lContainerPIDVars[2] = lRapOmega          ;
		lContainerPIDVars[3] = nTrackPrimaryMultiplicity;       // FIXME : nTrackPrimaryMultiplicity not set for AOD ... = -1
			
		// No PID
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorKaonForTPC  )
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
		if( lIsBachelorKaonForTPC && 
		    lIsPosProtonForTPC     )
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorKaonForTPC && 
		    lIsPosProtonForTPC    && 
		    lIsNegPionForTPC       )
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorKaon        )
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorKaon       && 
		    lIsPosInOmegaProton    )
			fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if(lIsBachelorKaon     && 
		   lIsPosInOmegaProton && 
		   lIsNegInOmegaPion    )
		 	fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Omega Plus		
	if( lChargeXi > 0 ) {
		lContainerPIDVars[0] = lXiTransvMom       ;
		lContainerPIDVars[1] = lInvMassOmegaPlus ;
		lContainerPIDVars[2] = lRapOmega          ;
		lContainerPIDVars[3] = nTrackPrimaryMultiplicity;       // FIXME : nTrackPrimaryMultiplicity not set for AOD ... = -1
			
		// No PID
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorKaonForTPC  )
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
		if( lIsBachelorKaonForTPC && 
		    lIsNegProtonForTPC     )
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorKaonForTPC && 
		    lIsNegProtonForTPC    && 
		    lIsPosPionForTPC       )
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorKaon        )
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorKaon       && 
		    lIsNegInOmegaProton    )
			fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if(lIsBachelorKaon     && 
		   lIsNegInOmegaProton && 
		   lIsPosInOmegaPion    )
		 	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
        	
        // II.Fill.Step 7 : filling the AliCFContainer (optimisation of topological selections)
        Double_t lContainerCutVars[20] = {0.0};
                        
        lContainerCutVars[0]  = lDcaXiDaughters;
        lContainerCutVars[1]  = lDcaBachToPrimVertexXi;
        lContainerCutVars[2]  = lXiCosineOfPointingAngle;
        lContainerCutVars[3]  = lXiRadius2D;
        lContainerCutVars[4]  = lInvMassLambdaAsCascDghter;
        lContainerCutVars[5]  = lDcaV0DaughtersXi;
        lContainerCutVars[6]  = lV0CosineOfPointingAngleXi;
        lContainerCutVars[7]  = lV0RadiusXi;
        lContainerCutVars[8]  = lDcaV0ToPrimVertexXi;	
        lContainerCutVars[9]  = lDcaPosToPrimVertexXi;
        lContainerCutVars[10] = lDcaNegToPrimVertexXi;
        
        lContainerCutVars[13] = lXiTransvMom;
        
        lContainerCutVars[16] = lBestPrimaryVtxPos[2];
        lContainerCutVars[17] = nTrackPrimaryMultiplicity;       // FIXME : nTrackPrimaryMultiplicity not checked for AOD ... = 0
        lContainerCutVars[18] = lSPDTrackletsMultiplicity;       // FIXME : SPDTrackletsMultiplicity is not available for AOD ... = -1
        lContainerCutVars[19] = lBachTPCClusters;                // FIXME : BachTPCClusters          is not available for AOD ... = -1
        
        if( lChargeXi < 0 ) {
                lContainerCutVars[11] = lInvMassXiMinus;
                lContainerCutVars[12] = 1.63;
                lContainerCutVars[14] = lRapXi;
                lContainerCutVars[15] = -1.;
                        if( lIsBachelorPionForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )       fCFContCascadeCuts->Fill(lContainerCutVars,0); // for Xi-

                lContainerCutVars[11] = 1.26;
                lContainerCutVars[12] = lInvMassOmegaMinus;
                lContainerCutVars[14] = -1.;
                lContainerCutVars[15] = lRapOmega;
                        if( lIsBachelorKaonForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )       fCFContCascadeCuts->Fill(lContainerCutVars,2); // for Omega-
        }
        else{
                lContainerCutVars[11] = lInvMassXiPlus;
                lContainerCutVars[12] = 1.63;
                lContainerCutVars[14] = lRapXi;
                lContainerCutVars[15] = -1.; 
                        if( lIsBachelorPionForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )       fCFContCascadeCuts->Fill(lContainerCutVars,1); // for Xi+

                lContainerCutVars[11] = 1.26;
                lContainerCutVars[12] = lInvMassOmegaPlus;
                lContainerCutVars[14] = -1.;
                lContainerCutVars[15] = lRapOmega;
                        if( lIsBachelorKaonForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )       fCFContCascadeCuts->Fill(lContainerCutVars,3); // for Omega+
        }

                        
        // II.Fill.Step 8 :  angular correlations
        
        if( lChargeXi < 0 ){
                if( lIsBachelorPionForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )      DoAngularCorrelation("Xi-",    lInvMassXiMinus,    lArrTrackID, lTVect3MomXi, lEta, lRapXi);
                if( lIsBachelorKaonForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )      DoAngularCorrelation("Omega-", lInvMassOmegaMinus, lArrTrackID, lTVect3MomXi, lEta, lRapOmega);
        }
        else{
                if( lIsBachelorPionForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )      DoAngularCorrelation("Xi+",    lInvMassXiPlus,    lArrTrackID, lTVect3MomXi, lEta, lRapXi);
                if( lIsBachelorKaonForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )      DoAngularCorrelation("Omega+", lInvMassOmegaPlus, lArrTrackID, lTVect3MomXi, lEta, lRapOmega);
        }
        
	
    }// end of the Cascade loop (ESD or AOD)
    
  
  // Post output data.
 PostData(1, fListHistCascade);
}


void AliAnalysisTaskCheckCascade::DoAngularCorrelation( const Char_t   *lCascType, 
							      Double_t  lInvMassCascade, 
							const Int_t    *lArrTrackID,
							      TVector3 &lTVect3MomXi, 
							      Double_t  lEtaXi,
                                                              Double_t  lRapCascade){
  // Perform the Delta(Phi)Delta(Eta) analysis 
  // by properly filling the THnSparseF 

        if( fAnalysisType == "AOD") return; // FIXME : AOD development lost, because of AliESDtrack needed by AliESDtrackCuts

        TString lStrCascType( lCascType );

        // Check the Xi- candidate is within the proper mass window m0 +- 10 MeV
        Double_t lCascPdgMass = 0.0;
        if( lStrCascType.Contains("Xi") )       lCascPdgMass = 1.3217;
        if( lStrCascType.Contains("Omega") )    lCascPdgMass = 1.6724;

        if( lInvMassCascade > lCascPdgMass + 0.010) return;
        if( lInvMassCascade < lCascPdgMass - 0.010) return;

        // Check the Xi- candidate is within the proper rapidity window (flat efficiency)
        if( lRapCascade >  0.5 ) return;
        if( lRapCascade < -0.5 ) return;


        // 3 options to follow for the correlations:
        // 1.1 - "TrigAnyCasc-AssoAnyPrim"
        // 1.2 - "TrigCascLeading-AssoAnyPrim"
        // 2.  - "TrigLeadingTrck-AssoCasc"

        if(fAngularCorrelationType.Contains("AssoAnyPrim") ){
        //----------------------- Option 1 ---------------------------------------------------------------------------------------------------
        // Cascade    = trigger, 
        // Associated = all the primary tracks in the event.
        //     1.1 Cascade = trigger above a certain pt but nothing more complicated than that
        //     1.2 Cascade = leading particle (à la heavy-ion -> test for coalescence)

                Bool_t kRejectLowPtCascades = kTRUE;
                if(kRejectLowPtCascades && (lTVect3MomXi.Pt() < 1.7) ) return;
                // Do not even consider the cascade of low pt for the correlation ...

                if(fAngularCorrelationType == "TrigCascLeading-AssoAnyPrim"){// Require the Cascade To be the Leading Part. in the event

                        // 1st loop: check there is no primary track with a higher pt ...
                        // = The cascade is meant to be a leading particle : Pt(Casc) > any primary track in the event
                        for(Int_t TrckIdx = 0; TrckIdx < (InputEvent())->GetNumberOfTracks() ; TrckIdx++ )
                        {// Loop over all the tracks of the event

                                AliESDtrack *lCurrentTrck = dynamic_cast<AliESDtrack*>( (InputEvent())->GetTrack( TrckIdx ) );
                                        if (!lCurrentTrck ) {
                                                AliWarning("ERROR Correl. Study : Could not retrieve a track while looping over the event tracks ...");
                                                continue;
                                        }

                                if( !fESDtrackCuts->AcceptTrack(lCurrentTrck) ) continue;
                                // Just consider primary tracks (= reject track that are not primary ones)
                                if(lTVect3MomXi.Pt() < lCurrentTrck->Pt() ) return;	
                                // Room for improvement: //FIXME
                                // 1. There is a given resolution on pt : maybe release the cut Pt(casc) < Pt(track)*90% ?
                                // 2. Apply this cut only when DeltaPhi(casc, track) > 90 deg = when track is in the near-side ?
                                // 3. Anti-splitting cut (like in Femto analysis) ? = now done via ESDtrackCuts ...

                        }// end control loop
                }// end of prelim. check : Cascade = leading part in the event

                // 2nd loop: filling loop
                for(Int_t TrckIdx = 0; TrckIdx < (InputEvent())->GetNumberOfTracks() ; TrckIdx++ )
                {// Loop over all the tracks of the event

                        AliESDtrack *lCurrentTrck = dynamic_cast<AliESDtrack*>( (InputEvent())->GetTrack( TrckIdx ) );
                                if (!lCurrentTrck ) {
                                        AliWarning("ERROR Correl. Study : Could not retrieve a track while looping over the event tracks ...");
                                        continue;
                                }
                        // Just consider primary tracks (= reject track that are not primary ones)
                        if( !fESDtrackCuts->AcceptTrack(lCurrentTrck) ) continue;

                        // Room for improvement: //FIXME
                        // 1.	
                        // 2. Exclude the tracks that build the condisdered cascade = the bachelor + the V0 dghters
                        //     This may bias the outcome, especially for low multplicity events.
                        // Note : For ESD event, track ID == track index.
                                if(lCurrentTrck->GetID() == lArrTrackID[0]) continue;
                                if(lCurrentTrck->GetID() == lArrTrackID[1]) continue;
                                if(lCurrentTrck->GetID() == lArrTrackID[2]) continue;

                        TVector3 lTVect3MomTrck(lCurrentTrck->Px(), lCurrentTrck->Py(), lCurrentTrck->Pz() );

                        // 2 hypotheses made here :
                        //   - The Xi trajectory is a straight line,
                        //   - The Xi doesn't loose any energy by crossing the first layer(s) of ITS, if ever;
                        //      So, meaning hyp: vect p(Xi) at the emission = vect p(Xi) at the decay vertex
                        //      By doing this, we introduce a systematic error on the cascade Phi ...
                        // Room for improvement: take into account the curvature of the Xi trajectory ?
                        //                       or rather, the resolution in space of the decay vertex ...
                        //FIXME

                        Double_t lHnSpFillVar[5] = {0.};
                        lHnSpFillVar[0] = lTVect3MomXi.DeltaPhi(lTVect3MomTrck) * 180.0/TMath::Pi(); // Delta phi(Casc,Track) (deg)
                        if(lHnSpFillVar[0] < -50.0) lHnSpFillVar[0] += 360.0;	
                        lHnSpFillVar[1] = lEtaXi - lCurrentTrck->Eta();                 // Delta eta(Casc,Track)
                        lHnSpFillVar[2] = lTVect3MomXi.Pt();                            // Pt_{Casc}
                        lHnSpFillVar[3] = lCurrentTrck->Pt();                           // Pt_{any track}
                        lHnSpFillVar[4] = lInvMassCascade;                              // Eff. Inv Mass (control var)

                        if(      lStrCascType.Contains("Xi-") )      fHnSpAngularCorrXiMinus    ->Fill( lHnSpFillVar );
                        else if( lStrCascType.Contains("Xi+") )      fHnSpAngularCorrXiPlus     ->Fill( lHnSpFillVar );
                        else if( lStrCascType.Contains("Omega-") )   fHnSpAngularCorrOmegaMinus ->Fill( lHnSpFillVar );
                        else if( lStrCascType.Contains("Omega+") )   fHnSpAngularCorrOmegaPlus  ->Fill( lHnSpFillVar );

                }// end - Loop over all the tracks in the event

        }// end of correlation type : "Trig = Casc - Asso = AnyPrim", (cases 1.1 and 1.2)



        else if(fAngularCorrelationType == "TrigLeadingTrck-AssoCasc"){

        //----------------------- Option 2 ---------------------------------------------------------------------------------------------------
        // Trigger    = trigger,
        // Associated = the cascade'S' in the event
        // NOTE : several good cascades could be present in the event (e.g. one leading track as trigger, 2 associated Xi) ...
        //       The present function will then be called several times.
        //        = issue for the normalisation ...

                // 1st loop: 
                //            find the index of the (1) primary track (2) which is the leading particle in pt
                // NOTE :  we do not take into account the Cascade pt, i.e. pt(Casc) could be greater or lower than pt(Leading) ...
                Int_t    lLeadingPrimTrackIdx = -1;
                Double_t lPtMax               =  0.1;

                for(Int_t TrckIdx = 0; TrckIdx < (InputEvent())->GetNumberOfTracks() ; TrckIdx++ )
                {// Loop over all the tracks of the event

                        AliESDtrack *lCurrentTrck = dynamic_cast<AliESDtrack*>( (InputEvent())->GetTrack( TrckIdx ) );
                                if (!lCurrentTrck ) {
                                        AliWarning("ERROR Correl. Study : Could not retrieve a track while looping over the event tracks ...");
                                        continue;
                                }
        
                        // Primary track selection
                        if( !fESDtrackCuts->AcceptTrack(lCurrentTrck) ) continue;
        
                        // Exclude the tracks that build the condisdered cascade = the bachelor + the V0 dghters
                        //     This may bias the outcome, especially for low multplicity events.
                        // Note : For ESD event, track ID == track index.
                        if(lCurrentTrck->GetID() == lArrTrackID[0]) continue;
                        if(lCurrentTrck->GetID() == lArrTrackID[1]) continue;
                        if(lCurrentTrck->GetID() == lArrTrackID[2]) continue;
        
                        // Towards the leading track
                        if( lPtMax < lCurrentTrck->Pt() ){
                                lLeadingPrimTrackIdx = TMath::Abs( lCurrentTrck->GetID() );
                                lPtMax               = lCurrentTrck->Pt();
                        }
                }// end leading track finding loop
        
                if( lLeadingPrimTrackIdx < 0 ) return;
                if( lPtMax < 0.101 )           return;
        
        
                // 2nd step: filling ONCE the THnSparse
                AliESDtrack *lLeadingTrck = dynamic_cast<AliESDtrack*>( (InputEvent())->GetTrack( lLeadingPrimTrackIdx ) );
                
                TVector3 lTVect3MomLeadingTrck( lLeadingTrck->Px(), lLeadingTrck->Py(), lLeadingTrck->Pz() );
                
                // 2 hypotheses made here :
                //   - The Xi trajectory is a straight line,
                //   - The Xi doesn't loose any energy by crossing the first layer(s) of ITS, if ever;
                //      So, meaning hyp: vect p(Xi) at the emission = vect p(Xi) at the decay vertex
                //      By doing this, we introduce a systematic error on the cascade Phi ...
                // Room for improvement: take into account the curvature of the Xi trajectory ?
                //                       or rather, the resolution in space of the decay vertex ...
                //FIXME
        
                Double_t lHnSpFillVar[5] = {0.};
                lHnSpFillVar[0] = lTVect3MomLeadingTrck.DeltaPhi(lTVect3MomXi) * 180.0/TMath::Pi(); // Delta phi(leading Track, Casc) (deg)
                if(lHnSpFillVar[0] < -50.0) lHnSpFillVar[0] += 360.0;
                lHnSpFillVar[1] = lLeadingTrck->Eta() - lEtaXi;                 // Delta eta(leading Track, Casc)
                lHnSpFillVar[2] = lTVect3MomXi.Pt();                            // Pt_{Casc}
                lHnSpFillVar[3] = lLeadingTrck->Pt();                           // Pt_{leading track}
                lHnSpFillVar[4] = lInvMassCascade;                              // Eff. Inv Mass (control var)
                
                if(      lStrCascType.Contains("Xi-") )      fHnSpAngularCorrXiMinus    ->Fill( lHnSpFillVar );
                else if( lStrCascType.Contains("Xi+") )      fHnSpAngularCorrXiPlus     ->Fill( lHnSpFillVar );
                else if( lStrCascType.Contains("Omega-") )   fHnSpAngularCorrOmegaMinus ->Fill( lHnSpFillVar );
                else if( lStrCascType.Contains("Omega+") )   fHnSpAngularCorrOmegaPlus  ->Fill( lHnSpFillVar );
        }// end of correlation type : "Trig = LeadingTrck -Asso = Casc"
        else
                return;

}

Int_t AliAnalysisTaskCheckCascade::DoESDTrackWithTPCrefitMultiplicity(const AliESDEvent *lESDevent)
{
    // Checking the number of tracks with TPCrefit for each event
    // Needed for a rough assessment of the event multiplicity
        
        Int_t nTrackWithTPCrefitMultiplicity = 0;
        for(Int_t iTrackIdx = 0; iTrackIdx < (InputEvent())->GetNumberOfTracks(); iTrackIdx++){
                AliESDtrack *esdTrack	= 0x0;
                             esdTrack	= lESDevent->GetTrack( iTrackIdx );
                if (!esdTrack) { AliWarning("Pb / Could not retrieve one track within the track loop for TPCrefit check ..."); continue; }

                ULong_t lTrackStatus    = esdTrack->GetStatus();
                            if ((lTrackStatus&AliESDtrack::kTPCrefit)    == 0) continue;
                            else nTrackWithTPCrefitMultiplicity++;
                    // FIXME :
                    // The goal here is to get a better assessment of the event multiplicity.
                    // (InputEvent())->GetNumberOfTracks() takes into account ITS std alone tracks + global tracks
                    // This may introduce a bias. Hence the number of TPC refit tracks.
                    // Note : the event multiplicity = analysis on its own... See Jacek's or Jan Fiete's analysis on dN/d(eta)

        }// end loop over all event tracks
        return  nTrackWithTPCrefitMultiplicity;
}





//________________________________________________________________________
void AliAnalysisTaskCheckCascade::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  TList *cRetrievedList = 0x0;
         cRetrievedList = (TList*)GetOutputData(1);
	if(!cRetrievedList){
		AliWarning("ERROR - AliAnalysisTaskCheckCascade: ouput data container list not available\n"); return;
	}

  fHistTrackMultiplicityForTrigEvt = dynamic_cast<TH1F*> (   cRetrievedList->FindObject("fHistTrackMultiplicityForTrigEvt") );
  if (!fHistTrackMultiplicityForTrigEvt) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascade: fHistTrackMultiplicityForTrigEvt not available\n"); return;
	}
 
  fHistCascadeMultiplicityForTrigEvt = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistCascadeMultiplicityForTrigEvt") );
	if (!fHistCascadeMultiplicityForTrigEvt) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascade: fHistCascadeMultiplicityForTrigEvt not available\n"); return;
	}
	
  fHistMassXiMinus    = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassXiMinus") );	
	if (!fHistMassXiMinus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascade: fHistMassXiMinus not available\n"); return;
	}
  fHistMassXiPlus     = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassXiPlus") );
	if (!fHistMassXiPlus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascade: fHistMassXiPlus not available\n"); return;
	}	
  fHistMassOmegaMinus = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassOmegaMinus") );
	if (!fHistMassOmegaMinus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascade: fHistMassOmegaMinus not available\n"); return;
	}
  fHistMassOmegaPlus  = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassOmegaPlus") );	
	if (!fHistMassOmegaPlus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascade: fHistMassOmegaPlus not available\n"); return;
	}
  
  TCanvas *canCheckCascade = new TCanvas("AliAnalysisTaskCheckCascade","CheckCascade overview",10,10,1010,660);
  canCheckCascade->Divide(2,2);
  
  canCheckCascade->cd(1);
  canCheckCascade->cd(1)->SetLogy();
  fHistTrackMultiplicityForTrigEvt->SetMarkerStyle(kFullStar);  
  fHistTrackMultiplicityForTrigEvt->GetXaxis()->SetLabelFont(42);
  fHistTrackMultiplicityForTrigEvt->GetYaxis()->SetLabelFont(42);
  fHistTrackMultiplicityForTrigEvt->SetTitleFont(42, "xy");
  fHistTrackMultiplicityForTrigEvt->GetXaxis()->SetTitleOffset(1.1);
  fHistTrackMultiplicityForTrigEvt->DrawCopy("H");
  
  canCheckCascade->cd(2);  
  canCheckCascade->cd(2)->SetLogy();
  fHistCascadeMultiplicityForTrigEvt->SetMarkerStyle(kOpenSquare);
  fHistCascadeMultiplicityForTrigEvt->GetXaxis()->SetLabelFont(42);
  fHistCascadeMultiplicityForTrigEvt->GetYaxis()->SetLabelFont(42);
  fHistCascadeMultiplicityForTrigEvt->SetTitleFont(42, "xy");
  fHistCascadeMultiplicityForTrigEvt->GetXaxis()->SetTitleOffset(1.1);
  fHistCascadeMultiplicityForTrigEvt->DrawCopy("E");
  
  canCheckCascade->cd(3);  
  fHistMassXiMinus ->SetMarkerStyle(kFullCircle);
  fHistMassXiMinus ->SetMarkerSize(0.5);
  fHistMassXiMinus ->GetXaxis()->SetLabelFont(42);
  fHistMassXiMinus ->GetYaxis()->SetLabelFont(42);
  fHistMassXiMinus ->SetTitleFont(42, "xy");
  fHistMassXiMinus ->GetXaxis()->SetTitleOffset(1.1);
  fHistMassXiMinus ->GetYaxis()->SetTitleOffset(1.3);
   // fHistMassXiMinus->Rebin(2);
  fHistMassXiMinus ->GetXaxis()->SetRangeUser(1.24, 1.42);
  fHistMassXiMinus ->DrawCopy("E");
  
  fHistMassXiPlus ->SetMarkerStyle(kOpenCircle);
  fHistMassXiPlus ->SetMarkerColor(kRed+2);
  fHistMassXiPlus ->SetLineColor(kRed+2);
  fHistMassXiPlus ->SetMarkerSize(0.5);
  // fHistMassXiPlus ->Rebin(2);
  fHistMassXiPlus ->DrawCopy("ESAME");
  
  
  TLegend *legendeXi =new TLegend(0.67,0.34,0.97,0.54);
 		legendeXi->SetTextFont(42);
 		legendeXi->SetTextSize(0.05);
 		legendeXi->SetFillColor(kWhite);
 		legendeXi->AddEntry( fHistMassXiMinus,"#Xi^{-} candidates","lp");
 		legendeXi->AddEntry( fHistMassXiPlus,"#Xi^{+} candidates","lp");
 		legendeXi->Draw();
  
  
  canCheckCascade->cd(4);  
  fHistMassOmegaPlus ->SetMarkerStyle(kOpenCircle);
  fHistMassOmegaPlus ->SetMarkerColor(kRed+2);
  fHistMassOmegaPlus ->SetLineColor(kRed+2);
  fHistMassOmegaPlus ->SetMarkerSize(0.5);
  fHistMassOmegaPlus ->GetXaxis()->SetLabelFont(42);
  fHistMassOmegaPlus ->GetYaxis()->SetLabelFont(42);
  fHistMassOmegaPlus ->SetTitleFont(42, "xy");
  fHistMassOmegaPlus ->GetXaxis()->SetTitleOffset(1.1);
  fHistMassOmegaPlus ->GetYaxis()->SetTitleOffset(1.25);
  // fHistMassOmegaPlus ->Rebin(2);
  fHistMassOmegaPlus ->GetXaxis()->SetRangeUser(1.6, 1.84);
  fHistMassOmegaPlus ->DrawCopy("E");
  
  fHistMassOmegaMinus->SetMarkerStyle(kFullCircle);
  fHistMassOmegaMinus->SetMarkerSize(0.5);
  // fHistMassOmegaMinus->Rebin(2);
  fHistMassOmegaMinus->DrawCopy("ESAME");

  
   TLegend *legendeOmega = new TLegend(0.67,0.34,0.97,0.54);
 		legendeOmega->SetTextFont(42);
 		legendeOmega->SetTextSize(0.05);
 		legendeOmega->SetFillColor(kWhite);
 		legendeOmega->AddEntry( fHistMassOmegaMinus,"#Omega^{-} candidates","lp");
 		legendeOmega->AddEntry( fHistMassOmegaPlus,"#Omega^{+} candidates","lp");
 		legendeOmega->Draw();

}
