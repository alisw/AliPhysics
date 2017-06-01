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
//            AliAnalysisTaskCheckCascadepp class
//
//            Origin AliAnalysisTaskCheckCascade which has four roles :
//              1. QAing the Cascades from ESD and AOD
//                 Origin:  AliAnalysisTaskESDCheckV0 by Boris Hippolyte Nov2007, hippolyt@in2p3.fr
//              2. Prepare the plots which stand as raw material for yield extraction (wi/wo PID)
//              3. Supply an AliCFContainer meant to define the optimised topological selections
//              4. Rough azimuthal correlation study (Eta, Phi)
//              Adapted to Cascade : A.Maire Mar2008, antonin.maire@ires.in2p3.fr
//              Modified :           A.Maire Mar2010 
//
//              Adapted to PbPb analysis: M. Nicassio, maria.nicassio@ba.infn.it
//               Feb-August2011
//                - Physics selection moved to the run.C macro
//                - Centrality selection added (+ setters) and histos
//                - flag and setters added (CF container usage, vertex range)
//                - histo added and histo/container binning changed 
//                - protection in the destructor for CAF usage          
//                - AliWarning disabled
//                - number of tracklets from AOD also          
//                - automatic settings for PID
//               September2011
//                - proper time histos/container added (V0 and Cascades)
//                - cosine PA V0 wrt Xi vertex in the container  
//               November2011
//                - re-run V0's and cascade's vertexers (SetCuts instead SetDefaultCuts!!)
//                - problems of libraries on Grid --> code copied in the task (from AliRoot v5-10-AN
//                  where new pt dependent V0's cosPA cut implemented by Iouri) 
//                - AOD analysis part completed 
//
//              Adapted to pp 2.76 TeV analysis: D. Colella, domenico.colella@ba.infn.it
//               Gen-now 2012
//                - Physics selection re-moved here (mainly for normalization in the efficiency calcuation)
//                - Centrality selection deleted
//
//              Adapted to pPb 5.02 TeV analysis: D. Colella, domenico.colella@ba.infn.it
//               Aug-Sep 2014
//               - Added the parameter fCollidingSystem, to distingish between pp and pPb procedures
//               Aug 2015
//               - Clarify the usage of SDD selection after the introduction of fCollifingSystem
//
//              Improvement for pp 13 TeV analysis: D. Colella, domenico.colella@ba.infn.it
//               Feb 2015
//               - Added the possibility to select events in the class kINT7 for the pp@13TeV analysis through "fkSwitchINT7"
//               - Revision of the event selections for the pp@13TeV analysis
//               May 2017
//               - 
//             
//
//-----------------------------------------------------------------

class TTree;
class TParticle;
class TVector3;

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
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"

#include "AliESDVZERO.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h" 
#include "AliAODInputHandler.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliAODTrack.h"
#include "AliAnalysisUtils.h"


#include "AliAnalysisTaskCheckCascadepp.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskCheckCascadepp)


//________________________________________________________________________
AliAnalysisTaskCheckCascadepp::AliAnalysisTaskCheckCascadepp() 
  : AliAnalysisTaskSE(), 
    fAnalysisType                   ("ESD"),
    fCollidingSystem                (0),
    fkTriggerClass                  (AliVEvent::kINT7),
    fApplyEvSelSDDstatus            (kFALSE),
    fApplyEvSelDAQincomplete        (kTRUE), 
    fApplyEvSelSPDclustervstracklet (kTRUE),
    fApplyEvSelPileup               (kTRUE),
    fApplyEvSelPhysicsSel           (kTRUE),
    fApplyEvSelNoTPConlyPrimVtx     (kTRUE),
    fApplyEvSelSPDvtxres            (kTRUE),
    fApplyEvSelVtxProximity         (kTRUE),
    fApplyEvSelZprimVtxPos          (kTRUE),
    fESDtrackCuts                   (0), 
    fUtils                          (0),
    fPIDResponse                    (0),
    fRerunV0CascVertexers           (kFALSE),
    fwithSDD                        (kFALSE),
    fExtraSelections                (kFALSE),
    fTrackQualityCutTPCrefit        (kTRUE),
    fTrackQualityCutnTPCcls         (kTRUE),
    fMinnTPCcls                     (0),
    fMinTPCcrossrawoverfindable     (0),
    fVtxRangeMax                    (10),
    fVtxRangeMin                    (0),
    fMinPtCutOnDaughterTracks       (0),
    fEtaCutOnDaughterTracks         (0),
    fSPDPileUpminContributors       (3),
    fTPCPIDsigma                    (4),
    fSuffix                         (""),

    // - Plots initialisation
    fListHistCascade(0),
      // Cascades multiplicity plots
      fHistCascadeMultiplicityBeforeAnySel(0),
      fHistCascadeMultiplicityAfterSDDstatusSel(0),
      fHistCascadeMultiplicityAfterDAQincompleteEvRej(0),
      fHistCascadeMultiplicityAfterSPDclustervstrackletSel(0),
      fHistCascadeMultiplicityAfterPileupRej(0),
      fHistCascadeMultiplicityAfterPhysicsSel(0),
      fHistCascadeMultiplicityAfterRevertexing(0),
      fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel(0),      
      fHistCascadeMultiplicityAfterSPDresolution(0),
      fHistCascadeMultiplicityAfterVerticesProximity(0),
      fHistCascadeMultiplicityAfterZprimVtxPosSel(0),           
      // Vertex position plots (BestVertex)
      fHistPVx(0), fHistPVy(0), fHistPVz(0),
      fHistPVxAnalysis(0), fHistPVyAnalysis(0), fHistPVzAnalysis(0),    
      // TPC cluster distributions for daughters
      fHistPosV0TPCClusters(0), 
      fHistNegV0TPCClusters(0), 
      fHistBachTPCClusters(0),
      // Cut's variables distributions
      fHistEffMassXi(0), 
      fHistDcaXiDaughters(0), 
      fHistDcaBachToPrimVertex(0), 
      fHistXiCosineOfPointingAngle(0), 
      fHistXiRadius(0),
      fHistMassLambdaAsCascDghter(0),
      fHistDcaV0DaughtersXi(0),
      fHistDcaV0ToPrimVertexXi(0), 
      fHistV0CosineOfPointingAngleXi(0),
      fHistV0RadiusXi(0),
      fHistDcaPosToPrimVertexXi(0), 
      fHistDcaNegToPrimVertexXi(0), 
      // Invariant mass distributions
      fHistMassXiMinus(0), fHistMassXiPlus(0), fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
      // Transverse and total momentum distributions
      fHistXiTransvMom(0), fHistXiTotMom(0), fHistBachTransvMomXi(0), fHistBachTotMomXi(0),
      // Others QA plots
      fHistChargeXi(0),
      fHistV0toXiCosineOfPointingAngle(0),
      fHistRapXi(0), fHistRapOmega(0), 
      fHistEtaXi(0), fHistEtaBachXi(0), fHistEtaPosXi(0), fHistEtaNegXi(0),
      fHistThetaXi(0), 
      fHistPhiXi(0),
      f2dHistArmenteros(0),			
      f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
      f2dHistEffMassLambdaVsEffMassXiPlus(0),  f2dHistEffMassXiVsEffMassOmegaPlus(0),
      f2dHistXiRadiusVsEffMassXiMinus(0),      f2dHistXiRadiusVsEffMassXiPlus(0),
      f2dHistXiRadiusVsEffMassOmegaMinus(0),   f2dHistXiRadiusVsEffMassOmegaPlus(0),
      f2dHistTPCdEdxOfCascDghters(0),
      f2dHistDcaXiDaughtersvsInvMass(0), 
      f2dHistDcaBachToPrimVertexvsInvMass(0), 
      f2dHistXiCosineOfPointingAnglevsInvMass(0),
      f2dHistMassLambdaAsCascDghtervsInvMass(0),
      f2dHistDcaV0DaughtersXivsInvMass(0),
      f2dHistDcaV0ToPrimVertexXivsInvMass(0),
      // Containers for cuts study 
      fCFContCascadePIDXiMinus(0),
      fCFContCascadePIDXiPlus(0),
      fCFContCascadePIDOmegaMinus(0),
      fCFContCascadePIDOmegaPlus(0),
      fCFContCascadeCuts(0)
    
    {
     // Dummy Constructor
        for(Int_t iV0selIdx   = 0; iV0selIdx   < 7; iV0selIdx++   ) { fV0Sels          [iV0selIdx   ] = -1.; }
        for(Int_t iCascSelIdx = 0; iCascSelIdx < 8; iCascSelIdx++ ) { fCascSels        [iCascSelIdx ] = -1.; }
    }


//________________________________________________________________________
AliAnalysisTaskCheckCascadepp::AliAnalysisTaskCheckCascadepp(const char *name) 
  : AliAnalysisTaskSE(name), 
    fAnalysisType                   ("ESD"),
    fCollidingSystem                (0),
    fkTriggerClass                  (AliVEvent::kINT7),
    fApplyEvSelSDDstatus            (kFALSE),
    fApplyEvSelDAQincomplete        (kTRUE),
    fApplyEvSelSPDclustervstracklet (kTRUE),
    fApplyEvSelPileup               (kTRUE),
    fApplyEvSelPhysicsSel           (kTRUE),
    fApplyEvSelNoTPConlyPrimVtx     (kTRUE),
    fApplyEvSelSPDvtxres            (kTRUE),
    fApplyEvSelVtxProximity         (kTRUE),
    fApplyEvSelZprimVtxPos          (kTRUE),
    fESDtrackCuts                   (0),
    fUtils                          (0),
    fPIDResponse                    (0),
    fRerunV0CascVertexers           (kFALSE),
    fwithSDD                        (kFALSE),
    fExtraSelections                (kFALSE),
    fTrackQualityCutTPCrefit        (kTRUE),
    fTrackQualityCutnTPCcls         (kTRUE),
    fMinnTPCcls                     (0),
    fMinTPCcrossrawoverfindable     (0),
    fVtxRangeMax                    (10),
    fVtxRangeMin                    (0),
    fMinPtCutOnDaughterTracks       (0),
    fEtaCutOnDaughterTracks         (0),
    fSPDPileUpminContributors       (3),
    fTPCPIDsigma                    (4),
    fSuffix                         (""),
   
    // - Plots initialisation
    fListHistCascade(0),
      // Cascades multiplicity plots
      fHistCascadeMultiplicityBeforeAnySel(0),
      fHistCascadeMultiplicityAfterSDDstatusSel(0),
      fHistCascadeMultiplicityAfterDAQincompleteEvRej(0),
      fHistCascadeMultiplicityAfterSPDclustervstrackletSel(0),
      fHistCascadeMultiplicityAfterPileupRej(0),
      fHistCascadeMultiplicityAfterPhysicsSel(0),
      fHistCascadeMultiplicityAfterRevertexing(0),
      fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel(0),   
      fHistCascadeMultiplicityAfterSPDresolution(0),
      fHistCascadeMultiplicityAfterVerticesProximity(0),
      fHistCascadeMultiplicityAfterZprimVtxPosSel(0),        
      // Vertex position plots (BestVertex)
      fHistPVx(0), fHistPVy(0), fHistPVz(0),
      fHistPVxAnalysis(0), fHistPVyAnalysis(0), fHistPVzAnalysis(0),
      // TPC cluster distributions for daughters
      fHistPosV0TPCClusters(0), fHistNegV0TPCClusters(0), fHistBachTPCClusters(0),
      // Cut's variables distributions
      fHistEffMassXi(0),
      fHistDcaXiDaughters(0),
      fHistDcaBachToPrimVertex(0),
      fHistXiCosineOfPointingAngle(0),
      fHistXiRadius(0),
      fHistMassLambdaAsCascDghter(0),
      fHistDcaV0DaughtersXi(0),
      fHistDcaV0ToPrimVertexXi(0),
      fHistV0CosineOfPointingAngleXi(0),
      fHistV0RadiusXi(0),
      fHistDcaPosToPrimVertexXi(0),
      fHistDcaNegToPrimVertexXi(0),
      // Invariant mass distributions
      fHistMassXiMinus(0), fHistMassXiPlus(0), fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
      // Transverse and total momentum distributions
      fHistXiTransvMom(0), fHistXiTotMom(0), fHistBachTransvMomXi(0), fHistBachTotMomXi(0),
      // Others QA plots
      fHistChargeXi(0),
      fHistV0toXiCosineOfPointingAngle(0),
      fHistRapXi(0), fHistRapOmega(0),
      fHistEtaXi(0), fHistEtaBachXi(0), fHistEtaPosXi(0), fHistEtaNegXi(0),
      fHistThetaXi(0),
      fHistPhiXi(0),
      f2dHistArmenteros(0),
      f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
      f2dHistEffMassLambdaVsEffMassXiPlus(0),  f2dHistEffMassXiVsEffMassOmegaPlus(0),
      f2dHistXiRadiusVsEffMassXiMinus(0),      f2dHistXiRadiusVsEffMassXiPlus(0),
      f2dHistXiRadiusVsEffMassOmegaMinus(0),   f2dHistXiRadiusVsEffMassOmegaPlus(0),
      f2dHistTPCdEdxOfCascDghters(0),
      f2dHistDcaXiDaughtersvsInvMass(0),
      f2dHistDcaBachToPrimVertexvsInvMass(0),
      f2dHistXiCosineOfPointingAnglevsInvMass(0),
      f2dHistMassLambdaAsCascDghtervsInvMass(0),
      f2dHistDcaV0DaughtersXivsInvMass(0),
      f2dHistDcaV0ToPrimVertexXivsInvMass(0),
      // Containers for cuts study 
      fCFContCascadePIDXiMinus(0),
      fCFContCascadePIDXiPlus(0),
      fCFContCascadePIDOmegaMinus(0),
      fCFContCascadePIDOmegaPlus(0),
      fCFContCascadeCuts(0)
    
    //_____Costructor____
    {
     // Define input and output slots here
     // Input slot #0 works with a TChain
     // DefineInput(0, TChain::Class());
     // Output slot #1 writes into a TList container (cascade)
     if (fCollidingSystem == 0) {
         fV0Sels[0] =  33.  ;     // max allowed chi2
         fV0Sels[1] =   0.01;     // min allowed impact parameter for the 1st daughter 
         fV0Sels[2] =   0.01;     // min allowed impact parameter for the 2nd daughter 
         fV0Sels[3] =   1.5;      // max allowed DCA between the daughter tracks       
         fV0Sels[4] =   0.9;      // min allowed cosine of V0's pointing angle - This is pT dependent        
         fV0Sels[5] =   0.2;      // min radius of the fiducial volume                 
         fV0Sels[6] = 200.;       // max radius of the fiducial volume                 
         fCascSels[0] =  33.;     // max allowed chi2 (same as PDC07)
         fCascSels[1] =   0.01;   // min allowed V0 impact parameter                    
         fCascSels[2] =   0.008;  // "window" around the Lambda mass                    
         fCascSels[3] =   0.01;   // min allowed bachelor's impact parameter          
         fCascSels[4] =   2.0;    // max allowed DCA between the V0 and the bachelor    
         fCascSels[5] =   0.95;   // min allowed cosine of the cascade pointing angle   
         fCascSels[6] =   0.2;    // min radius of the fiducial volume                  
         fCascSels[7] = 100.;     // max radius of the fiducial volume 
     } else if (fCollidingSystem == 1) {
         fV0Sels[0] =  33.  ;     // max allowed chi2
         fV0Sels[1] =   0.02;     // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
         fV0Sels[2] =   0.02;     // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
         fV0Sels[3] =   2.0 ;     // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
         fV0Sels[4] =   0.95;     // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
         fV0Sels[5] =   1.0 ;     // min radius of the fiducial volume                 (LHC09a4 : 0.2)
         fV0Sels[6] = 200.  ;     // max radius of the fiducial volume                 (LHC09a4 : 100.0)
         fCascSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
         fCascSels[1] =   0.05 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
         fCascSels[2] =   0.010;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
         fCascSels[3] =   0.03 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
         fCascSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
         fCascSels[5] =   0.95 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
         fCascSels[6] =   0.4  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
         fCascSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
     }
     // Output slot #0 writes into a TList container (Cascade)
     DefineOutput(1, TList::Class());
     DefineOutput(2, AliCFContainer::Class());
     DefineOutput(3, AliCFContainer::Class());
     DefineOutput(4, AliCFContainer::Class());
     DefineOutput(5, AliCFContainer::Class());
     DefineOutput(6, AliCFContainer::Class());
     AliLog::SetClassDebugLevel("AliAnalysisTaskCheckCascadepp",1);
    } 


    //_____Destructor_____
    AliAnalysisTaskCheckCascadepp::~AliAnalysisTaskCheckCascadepp() {
      // For all TH1, 2, 3 HnSparse and CFContainer are in the fListCascade TList.
      // They will be deleted when fListCascade is deleted by the TSelector dtor
      // Because of TList::SetOwner() ...   
       if (fListHistCascade && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())           { delete fListHistCascade; fListHistCascade = 0x0; }
       if (fCFContCascadePIDXiMinus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())   { delete fCFContCascadePIDXiMinus; fCFContCascadePIDXiMinus = 0x0; }
       if (fCFContCascadePIDXiPlus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())    { delete fCFContCascadePIDXiPlus; fCFContCascadePIDXiPlus = 0x0; }
       if (fCFContCascadePIDOmegaMinus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){ delete fCFContCascadePIDOmegaMinus; fCFContCascadePIDOmegaMinus = 0x0; }
       if (fCFContCascadePIDOmegaPlus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fCFContCascadePIDOmegaPlus; fCFContCascadePIDOmegaPlus = 0x0; }  
       if (fCFContCascadeCuts && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())         { delete fCFContCascadeCuts; fCFContCascadeCuts = 0x0; }
       if (fESDtrackCuts)                                                                          { delete fESDtrackCuts; fESDtrackCuts = 0x0; }
    }


//________________________________________________________________________
void AliAnalysisTaskCheckCascadepp::UserCreateOutputObjects() {
  // Create histograms
  // Called once
 fListHistCascade = new TList();
 fListHistCascade->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

 //-----------------------------------------------
 // Particle Identification Setup (new PID object)
 //-----------------------------------------------
 AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
 AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
 fPIDResponse = inputHandler->GetPIDResponse();
 // Only used to get the number of primary reconstructed tracks
 if (fAnalysisType == "ESD" && (! fESDtrackCuts )){
   fESDtrackCuts = new AliESDtrackCuts();
 }

 //---------------------------------------------------
 // Initialize cuts to re-run V0 and cascade vertexers
 //---------------------------------------------------
 // Not validated; to be checked
 if (fCollidingSystem == 0) {
      fV0Sels[0] =  33.  ;     // max allowed chi2
      fV0Sels[1] =   0.01;     // min allowed impact parameter for the 1st daughter 
      fV0Sels[2] =   0.01;     // min allowed impact parameter for the 2nd daughter 
      fV0Sels[3] =   1.5;      // max allowed DCA between the daughter tracks       
      fV0Sels[4] =   0.9;      // min allowed cosine of V0's pointing angle         
      fV0Sels[5] =   0.2;      // min radius of the fiducial volume                 
      fV0Sels[6] = 200.;       // max radius of the fiducial volume                 
      fCascSels[0] =  33.;     // max allowed chi2 (same as PDC07)
      fCascSels[1] =   0.01;   // min allowed V0 impact parameter                    
      fCascSels[2] =   0.008;  // "window" around the Lambda mass                    
      fCascSels[3] =   0.01;   // min allowed bachelor's impact parameter          
      fCascSels[4] =   2.0;    // max allowed DCA between the V0 and the bachelor    
      fCascSels[5] =   0.95;   // min allowed cosine of the cascade pointing angle   
      fCascSels[6] =   0.2;    // min radius of the fiducial volume                  
      fCascSels[7] = 100.;     // max radius of the fiducial volume 
 } else if (fCollidingSystem == 1) {
      fV0Sels[0] =  33.  ;     // max allowed chi2
      fV0Sels[1] =   0.02;     // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
      fV0Sels[2] =   0.02;     // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
      fV0Sels[3] =   2.0 ;     // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
      fV0Sels[4] =   0.95;     // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
      fV0Sels[5] =   1.0 ;     // min radius of the fiducial volume                 (LHC09a4 : 0.2)
      fV0Sels[6] = 200.  ;     // max radius of the fiducial volume                 (LHC09a4 : 100.0)
      fCascSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
      fCascSels[1] =   0.05 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
      fCascSels[2] =   0.010;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
      fCascSels[3] =   0.03 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
      fCascSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
      fCascSels[5] =   0.95 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
      fCascSels[6] =   0.4  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
      fCascSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
 }

 //----------------------
 // Initialize the histos
 //----------------------
 
 // - Cascades multiplicity plots 
 if(! fHistCascadeMultiplicityBeforeAnySel) {
        fHistCascadeMultiplicityBeforeAnySel = new TH1F("fHistCascadeMultiplicityBeforeAnySel",
                 "Cascades per event (before any selections);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityBeforeAnySel);
 }
 if(! fHistCascadeMultiplicityAfterSDDstatusSel) {
        fHistCascadeMultiplicityAfterSDDstatusSel = new TH1F("fHistCascadeMultiplicityAfterSDDstatusSel",
                 "Cascades per event (after the SDD status selection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterSDDstatusSel);
 }
 if(! fHistCascadeMultiplicityAfterDAQincompleteEvRej) {
        fHistCascadeMultiplicityAfterDAQincompleteEvRej = new TH1F("fHistCascadeMultiplicityAfterDAQincompleteEvRej",
                 "Cascades per event (after DAQ incomplete event rejection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterDAQincompleteEvRej);
 }
 if(! fHistCascadeMultiplicityAfterSPDclustervstrackletSel) {
        fHistCascadeMultiplicityAfterSPDclustervstrackletSel = new TH1F("fHistCascadeMultiplicityAfterSPDclustervstrackletSel",
                 "Cascades per event (after background rejection based on SPD cluster vs tracklet);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterSPDclustervstrackletSel);
 }
 if(! fHistCascadeMultiplicityAfterPileupRej) {
        fHistCascadeMultiplicityAfterPileupRej = new TH1F("fHistCascadeMultiplicityAfterPileupRej",
                 "Cascades per event (after pile-up events rejection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterPileupRej);
 }
 if(! fHistCascadeMultiplicityAfterPhysicsSel) {
        fHistCascadeMultiplicityAfterPhysicsSel = new TH1F("fHistCascadeMultiplicityAfterPhysicsSel",
                 "Cascades per event (after physics selection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterPhysicsSel);
 }
 if(! fHistCascadeMultiplicityAfterRevertexing) {
        fHistCascadeMultiplicityAfterRevertexing = new TH1F("fHistCascadeMultiplicityAfterRevertexing",
                 "Cascades per event (after revertexing);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterRevertexing);
 }
 if(! fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel) {
        fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel = new TH1F("fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel",
                 "Cascades per event (for selected events with well-established PV);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel);
 }
 if(! fHistCascadeMultiplicityAfterSPDresolution) {
        fHistCascadeMultiplicityAfterSPDresolution = new TH1F("fHistCascadeMultiplicityAfterSPDresolution",
                 "Cascades per event (after SPD resolution);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterSPDresolution);
 }
 if(! fHistCascadeMultiplicityAfterVerticesProximity) {
        fHistCascadeMultiplicityAfterVerticesProximity = new TH1F("fHistCascadeMultiplicityAfterVerticesProximity",
                 "Cascades per event (after vertices proximity);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterVerticesProximity);
 }
 if(! fHistCascadeMultiplicityAfterZprimVtxPosSel) {
        fHistCascadeMultiplicityAfterZprimVtxPosSel = new TH1F("fHistCascadeMultiplicityAfterZprimVtxPosSel",
                 "Cascades per event (after vertex cut selection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterZprimVtxPosSel);
 }
 // - Vertex position plots
 if(! fHistPVx ){
	fHistPVx = new TH1F("fHistPVx", "Best PV position in x; x (cm); Events", 2000, -0.5, 0.5);
	fListHistCascade->Add(fHistPVx);
 }
 if(! fHistPVy ){
        fHistPVy = new TH1F("fHistPVy", "Best PV position in y; y (cm); Events", 4000, -1.0, 1.0);
        fListHistCascade->Add(fHistPVy);
 }
 if(! fHistPVz ){
        fHistPVz = new TH1F("fHistPVz", "Best PV position in z; z (cm); Events", 400, -20, 20);
        fListHistCascade->Add(fHistPVz);
 } 
 if(! fHistPVxAnalysis ){
        fHistPVxAnalysis = new TH1F("fHistPVxAnalysis", "Best PV position in x (after events selections); x (cm); Events", 2000, -0.5, 0.5);
        fListHistCascade->Add(fHistPVxAnalysis);
 }
 if(! fHistPVyAnalysis ){
        fHistPVyAnalysis = new TH1F("fHistPVyAnalysis", "Best PV position in y (after events selections); y (cm); Events", 4000, -1.0, 1.0);
        fListHistCascade->Add(fHistPVyAnalysis);
 }
 if(! fHistPVzAnalysis ){
        fHistPVzAnalysis = new TH1F("fHistPVzAnalysis", "Best PV position in z (after events selections); z (cm); Events", 400, -20, 20);
        fListHistCascade->Add(fHistPVzAnalysis);
 }
 // - TPC clusetr sdistributions for daughters (histos for events containing at least ONE CASCADE)
 if(! fHistPosV0TPCClusters ){
        fHistPosV0TPCClusters = new TH1F("fHistPosV0TPCClusters", "TPC clusters for Pos. V0 daughter track, in Casc; Nbr of TPC clusters (V0 Pos.); Track counts", 165, 0.0, 165.0);
        fListHistCascade->Add(fHistPosV0TPCClusters);
 }
 if(! fHistNegV0TPCClusters ){
        fHistNegV0TPCClusters = new TH1F("fHistNegV0TPCClusters", "TPC clusters for Neg. V0 daughter track, in Casc; Nbr of TPC clusters (V0 Neg.); Track counts", 165, 0.0, 165.0);
        fListHistCascade->Add(fHistNegV0TPCClusters);
 }
 if(! fHistBachTPCClusters ){
        fHistBachTPCClusters = new TH1F("fHistBachTPCClusters", "TPC clusters for Bachelor track; Nbr of TPC clusters (Bach); Track counts", 165, 0.0, 165.0);
        fListHistCascade->Add(fHistBachTPCClusters);
 }
 // - Cut's variables distributions (typical histos for cascades): as example only for the Xi (both particle and anti-particle)
 if(! fHistEffMassXi) {
     fHistEffMassXi = new TH1F("fHistEffMassXi", "Xi candidates; Invariant Mass (GeV/c^{2}); Counts", 400, 1.2, 2.0);
     fListHistCascade->Add(fHistEffMassXi);
 }  
 if(! fHistDcaXiDaughters ){
     fHistDcaXiDaughters = new TH1F("fHistDcaXiDaughters", "DCA between Xi daughters; DCA (cm); Counts", 210, 0., 2.1);
     fListHistCascade->Add(fHistDcaXiDaughters);
 }
 if(! fHistDcaBachToPrimVertex) {
     fHistDcaBachToPrimVertex = new TH1F("fHistDcaBachToPrimVertex", "Impact parameter of Bach. to Prim. Vertex; DCA (cm); Counts", 250, 0., 0.25);
     fListHistCascade->Add(fHistDcaBachToPrimVertex);
 }
 if(! fHistXiCosineOfPointingAngle) {
     fHistXiCosineOfPointingAngle = new TH1F("fHistXiCosineOfPointingAngle", "Cosine of Xi Pointing Angle; Cos (Xi Point.Angl); Counts", 601, 0.94, 1.0001);
     fListHistCascade->Add(fHistXiCosineOfPointingAngle);
 }
 if(! fHistXiRadius ){
     fHistXiRadius = new TH1F("fHistXiRadius", "Cascade decay transv. radius; r (cm); Counts" , 2050, 0., 205.0);
     fListHistCascade->Add(fHistXiRadius);
 }
 if(! fHistMassLambdaAsCascDghter) {
     fHistMassLambdaAsCascDghter = new TH1F("fHistMassLambdaAsCascDghter", "#Lambda associated to cascade candidates; Eff. Mass (GeV/c^{2}); Counts", 300, 1.0, 1.3);
     fListHistCascade->Add(fHistMassLambdaAsCascDghter);
 }
 if(! fHistDcaV0DaughtersXi) {
     fHistDcaV0DaughtersXi = new TH1F("fHistDcaV0DaughtersXi", "DCA between V0 daughters, in cascade; DCA (cm); Counts", 320, 0., 1.6);
     fListHistCascade->Add(fHistDcaV0DaughtersXi);
 }
 if(! fHistDcaV0ToPrimVertexXi) {
     fHistDcaV0ToPrimVertexXi = new TH1F("fHistDcaV0ToPrimVertexXi", "Impact parameter of V0  to Prim. Vertex, in cascade; DCA (cm); Counts", 200, 0., 1.);
     fListHistCascade->Add(fHistDcaV0ToPrimVertexXi);
 }
 if(! fHistV0CosineOfPointingAngleXi) {
     fHistV0CosineOfPointingAngleXi = new TH1F("fHistV0CosineOfPointingAngleXi", "Cosine of V0 Pointing Angle, in cascade; Cos(V0 Point. Angl); Counts", 201, 0.8, 1.001);
     fListHistCascade->Add(fHistV0CosineOfPointingAngleXi);
 }
 if(! fHistV0RadiusXi) {
     fHistV0RadiusXi = new TH1F("fHistV0RadiusXi", "V0 decay radius, in cascade; radius (cm); Counts", 2050, 0., 205.0);
     fListHistCascade->Add(fHistV0RadiusXi);
 }
 if(! fHistDcaPosToPrimVertexXi) {
     fHistDcaPosToPrimVertexXi = new TH1F("fHistDcaPosToPrimVertexXi", "Impact parameter of V0 pos daughter to Prim. Vertex; DCA (cm); Counts", 300, 0, 3);
     fListHistCascade->Add(fHistDcaPosToPrimVertexXi);
 }
 if(! fHistDcaNegToPrimVertexXi) {
     fHistDcaNegToPrimVertexXi = new TH1F("fHistDcaNegToPrimVertexXi", "Impact parameter of V0 neg daughter to Prim. Vertex; DCA (cm); Counts", 300, 0, 3);
     fListHistCascade->Add(fHistDcaNegToPrimVertexXi);
 }
 // - Effective mass histos for cascades.
    //By cascade hyp  
 if(! fHistMassXiMinus) {
     fHistMassXiMinus = new TH1F("fHistMassXiMinus", "#Xi^{-} candidates; M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 400, 1.2, 2.0);
     fListHistCascade->Add(fHistMassXiMinus);
 } 
 if(! fHistMassXiPlus) {
     fHistMassXiPlus = new TH1F("fHistMassXiPlus", "#Xi^{+} candidates; M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts", 400, 1.2, 2.0);
     fListHistCascade->Add(fHistMassXiPlus);
 }
 if(! fHistMassOmegaMinus) {
     fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus", "#Omega^{-} candidates; M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500, 1.5, 2.5);
     fListHistCascade->Add(fHistMassOmegaMinus);
 }
 if(! fHistMassOmegaPlus) {
     fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus", "#Omega^{+} candidates;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2}); Counts", 500, 1.5, 2.5);
     fListHistCascade->Add(fHistMassOmegaPlus);
 }
 // - Transverse and total momentum distributions
 if(! fHistXiTransvMom ){
     fHistXiTransvMom = new TH1F("fHistXiTransvMom", "#Xi transverse momentum (cand. around the mass peak); p_{t}(#Xi) (GeV/c); Counts", 100, 0.0, 10.0);
     fListHistCascade->Add(fHistXiTransvMom);
 }
 if(! fHistXiTotMom ){
     fHistXiTotMom = new TH1F("fHistXiTotMom", "#Xi momentum norm (cand. around the mass peak); p_{tot}(#Xi) (GeV/c); Counts", 150, 0.0, 15.0);
     fListHistCascade->Add(fHistXiTotMom);
 }
 if(! fHistBachTransvMomXi ){
     fHistBachTransvMomXi = new TH1F("fHistBachTransvMomXi", "#Xi Bach. transverse momentum (cand. around the mass peak); p_{t}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
     fListHistCascade->Add(fHistBachTransvMomXi);
 }
 if(! fHistBachTotMomXi ){
     fHistBachTotMomXi = new TH1F("fHistBachTotMomXi", "#Xi Bach. momentum norm (cand. around the mass peak); p_{tot}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
     fListHistCascade->Add(fHistBachTotMomXi);
 }
 // - Others QA plots
    //TH1
 if(! fHistChargeXi ){
     fHistChargeXi = new TH1F("fHistChargeXi", "Charge of Xi candidates; Sign; Counts", 5, -2.0, 3.0);
     fListHistCascade->Add(fHistChargeXi);
 }
 if(! fHistV0toXiCosineOfPointingAngle) {
     fHistV0toXiCosineOfPointingAngle = new TH1F("fHistV0toXiCosineOfPointingAngle", "Cos. of V0 Ptng Angl / Xi vtx ; Cos(V0 Point. Angl / Xi vtx); Counts", 1101, 0.89, 1.0001);
     fListHistCascade->Add(fHistV0toXiCosineOfPointingAngle);
 }
 if(! fHistRapXi ){
     fHistRapXi = new TH1F("fHistRapXi", "Rapidity of #Xi candidates (around the mass peak); y; Counts", 20, -1.0, 1.0);
     fListHistCascade->Add(fHistRapXi);
 }
 if(! fHistRapOmega ){
     fHistRapOmega = new TH1F("fHistRapOmega", "Rapidity of #Omega candidates (around the mass peak); y; Counts", 20, -1.0, 1.0);
     fListHistCascade->Add(fHistRapOmega);
 }
 if(! fHistEtaXi ){
     fHistEtaXi = new TH1F("fHistEtaXi", "Pseudo-rap. of #Xi candidates (around the mass peak); #eta; Counts", 20, -1.0, 1.0);
     fListHistCascade->Add(fHistEtaXi);
 }
 if(! fHistEtaBachXi){
     fHistEtaBachXi = new TH1F("fHistEtaBachXi", "Pseudo-rap. of #Xi bachelor; #eta; Counts", 40, -2.0, 2.0);
     fListHistCascade->Add(fHistEtaBachXi);
 }
 if(! fHistEtaPosXi){
     fHistEtaPosXi = new TH1F("fHistEtaPosXi", "Pseudo-rap. of #Xi positive meson daughter; #eta; Counts", 40, -2.0, 2.0);
     fListHistCascade->Add(fHistEtaPosXi);
 }
 if(! fHistEtaNegXi){
     fHistEtaNegXi = new TH1F("fHistEtaNegXi", "Pseudo-rap. of #Xi negative meson daughter; #eta; Counts", 40, -2.0, 2.0);
     fListHistCascade->Add(fHistEtaNegXi);
 }
 if(! fHistThetaXi ){
     fHistThetaXi = new TH1F("fHistThetaXi", "#theta of #Xi candidates (around the mass peak); #theta (deg); Counts", 180, 0., 180.0);
     fListHistCascade->Add(fHistThetaXi);
 }
 if(! fHistPhiXi ){
     fHistPhiXi = new TH1F("fHistPhiXi", "#phi of #Xi candidates (around the mass peak); #phi (deg); Counts", 360, 0., 360.);
     fListHistCascade->Add(fHistPhiXi);
 }
 if(! f2dHistArmenteros) {
     f2dHistArmenteros = new TH2F("f2dHistArmenteros", "#alpha_{Arm}(casc. cand.) Vs Pt_{Arm}(casc. cand.); #alpha_{Arm}; Pt_{Arm} (GeV/c)", 140, -1.2, 1.2, 300, 0., 0.3);
     fListHistCascade->Add(f2dHistArmenteros);
 }
    //TH2
 if(! f2dHistEffMassLambdaVsEffMassXiMinus) {
     f2dHistEffMassLambdaVsEffMassXiMinus = new TH2F("f2dHistEffMassLambdaVsEffMassXiMinus", "M_{#Lambda} Vs M_{#Xi^{-} candidates}; Inv. M_{#Lambda^{0}} (GeV/c^{2}); M( #Lambda , #pi^{-} ) (GeV/c^{2})", 300, 1.1, 1.13, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistEffMassLambdaVsEffMassXiMinus);
 }
 if(! f2dHistEffMassXiVsEffMassOmegaMinus) {
     f2dHistEffMassXiVsEffMassOmegaMinus = new TH2F("f2dHistEffMassXiVsEffMassOmegaMinus", "M_{#Xi^{-} candidates} Vs M_{#Omega^{-} candidates}; M( #Lambda , #pi^{-} ) (GeV/c^{2}); M( #Lambda , K^{-} ) (GeV/c^{2})", 400, 1.2, 2.0, 500, 1.5, 2.5);
     fListHistCascade->Add(f2dHistEffMassXiVsEffMassOmegaMinus);
 }
 if(! f2dHistEffMassLambdaVsEffMassXiPlus) {
     f2dHistEffMassLambdaVsEffMassXiPlus = new TH2F("f2dHistEffMassLambdaVsEffMassXiPlus", "M_{#Lambda} Vs M_{#Xi^{+} candidates}; Inv. M_{#Lambda^{0}} (GeV/c^{2}); M( #Lambda , #pi^{+} ) (GeV/c^{2})", 300, 1.1, 1.13, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistEffMassLambdaVsEffMassXiPlus);
 }
 if(! f2dHistEffMassXiVsEffMassOmegaPlus) {
     f2dHistEffMassXiVsEffMassOmegaPlus = new TH2F("f2dHistEffMassXiVsEffMassOmegaPlus", "M_{#Xi^{+} candidates} Vs M_{#Omega^{+} candidates}; M( #Lambda , #pi^{+} ) (GeV/c^{2}); M( #Lambda , K^{+} ) (GeV/c^{2})", 400, 1.2, 2.0, 500, 1.5, 2.5);
     fListHistCascade->Add(f2dHistEffMassXiVsEffMassOmegaPlus);
 }
 if(! f2dHistXiRadiusVsEffMassXiMinus) {
     f2dHistXiRadiusVsEffMassXiMinus = new TH2F("f2dHistXiRadiusVsEffMassXiMinus", "Transv. R_{Xi Decay} Vs M_{#Xi^{-} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{-} ) (GeV/c^{2})", 450, 0., 45.0, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiMinus);
 }
 if(! f2dHistXiRadiusVsEffMassXiPlus) {
     f2dHistXiRadiusVsEffMassXiPlus = new TH2F("f2dHistXiRadiusVsEffMassXiPlus", "Transv. R_{Xi Decay} Vs M_{#Xi^{+} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{+} ) (GeV/c^{2})", 450, 0., 45.0, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiPlus);
 }
 if(! f2dHistXiRadiusVsEffMassOmegaMinus) {
     f2dHistXiRadiusVsEffMassOmegaMinus = new TH2F("f2dHistXiRadiusVsEffMassOmegaMinus", "Transv. R_{Xi Decay} Vs M_{#Omega^{-} candidates}; r_{cascade} (cm); M( #Lambda , K^{-} ) (GeV/c^{2}) ", 450, 0., 45.0, 500, 1.5, 2.5);
     fListHistCascade->Add(f2dHistXiRadiusVsEffMassOmegaMinus);
 }
 if(! f2dHistXiRadiusVsEffMassOmegaPlus) {
     f2dHistXiRadiusVsEffMassOmegaPlus = new TH2F("f2dHistXiRadiusVsEffMassOmegaPlus", "Transv. R_{Xi Decay} Vs M_{#Omega^{+} candidates}; r_{cascade} (cm); M( #Lambda , K^{+} ) (GeV/c^{2}) ", 450, 0., 45.0, 500, 1.5, 2.5);
     fListHistCascade->Add(f2dHistXiRadiusVsEffMassOmegaPlus);
 }
 if(! f2dHistTPCdEdxOfCascDghters){
     f2dHistTPCdEdxOfCascDghters = new TH2F("f2dHistTPCdEdxOfCascDghters", "TPC dE/dx of the cascade daughters; charge x || #vec{p}_{TPC inner wall}(Casc. daughter) || (GeV/c); TPC signal (ADC)", 2000, -10.0, 10.0, 450, 0., 900.);
     fListHistCascade->Add(f2dHistTPCdEdxOfCascDghters);
 }
 if(! f2dHistDcaXiDaughtersvsInvMass){
     f2dHistDcaXiDaughtersvsInvMass = new TH2F("f2dHistDcaXiDaughtersvsInvMass", "DCA between Xi Daughters; DCA (cm); Number of Cascades", 100, 0., 0.5, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistDcaXiDaughtersvsInvMass);
 }
 if(! f2dHistDcaBachToPrimVertexvsInvMass) {
     f2dHistDcaBachToPrimVertexvsInvMass = new TH2F("f2dHistDcaBachToPrimVertexvsInvMass", "DCA of Bach. to Prim. Vertex; DCA (cm); Number of Cascades", 250, 0., 0.25, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistDcaBachToPrimVertexvsInvMass);
 }
 if(! f2dHistXiCosineOfPointingAnglevsInvMass){
     f2dHistXiCosineOfPointingAnglevsInvMass = new TH2F("f2dHistXiCosineOfPointingAnglevsInvMass", "Cosine of Xi Pointing Angle; Cos (Xi Point.Angl); Number of Xis", 200, 0.99, 1.0, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistXiCosineOfPointingAnglevsInvMass);
 }
 if(! f2dHistMassLambdaAsCascDghtervsInvMass){ 
     f2dHistMassLambdaAsCascDghtervsInvMass = new TH2F("f2dHistMassLambdaAsCascDghtervsInvMass","#Lambda associated to Casc. candidates; Eff. Mass (GeV/c^{2}); Counts", 300, 1.00, 1.3, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistMassLambdaAsCascDghtervsInvMass);
 }
 if(! f2dHistDcaV0DaughtersXivsInvMass){
     f2dHistDcaV0DaughtersXivsInvMass = new TH2F("f2dHistDcaV0DaughtersXivsInvMass", "DCA between V0 daughters, in cascade; DCA (cm); Number of V0s", 120, 0., 0.6, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistDcaV0DaughtersXivsInvMass);
 }
 if(! f2dHistDcaV0ToPrimVertexXivsInvMass){
     f2dHistDcaV0ToPrimVertexXivsInvMass = new TH2F("f2dHistDcaV0ToPrimVertexXivsInvMass", "DCA of V0 to Prim. Vertex, in cascade; DCA (cm); Number of Cascades", 200, 0., 1., 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistDcaV0ToPrimVertexXivsInvMass);
 }
 // - Usefull string
 TString sddstatus = "";
 if      (fCollidingSystem == 0 && fApplyEvSelSDDstatus && fwithSDD)  sddstatus = "_wSDDon";
 else if (fCollidingSystem == 0 && fApplyEvSelSDDstatus && !fwithSDD) sddstatus = "_wSDDoff";
 TString cfcontname_cascpidximinus = Form("fCFContCascadePIDXiMinus_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fMinTPCcrossrawoverfindable,fVtxRangeMax,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks);
 cfcontname_cascpidximinus.Append(Form("%s",sddstatus.Data()));
 cfcontname_cascpidximinus.Append(Form("%s",fSuffix.Data()));
 TString cfcontname_cascpidxiplus = Form("fCFContCascadePIDXiPlus_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fMinTPCcrossrawoverfindable,fVtxRangeMax,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks);
 cfcontname_cascpidxiplus.Append(Form("%s",sddstatus.Data()));
 cfcontname_cascpidxiplus.Append(Form("%s",fSuffix.Data()));
 TString cfcontname_cascpidomegaminus = Form("fCFContCascadePIDOmegaMinus_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fMinTPCcrossrawoverfindable,fVtxRangeMax,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks);
 cfcontname_cascpidomegaminus.Append(Form("%s",sddstatus.Data()));
 cfcontname_cascpidomegaminus.Append(Form("%s",fSuffix.Data()));
 TString cfcontname_cascpidomegaplus = Form("fCFContCascadePIDOmegaPlus_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fMinTPCcrossrawoverfindable,fVtxRangeMax,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks);
 cfcontname_cascpidomegaplus.Append(Form("%s",sddstatus.Data()));
 cfcontname_cascpidomegaplus.Append(Form("%s",fSuffix.Data()));
 TString cfcontname_casccuts = Form("fCFContCascadeCuts_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fMinTPCcrossrawoverfindable,fVtxRangeMax,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks);
 cfcontname_casccuts.Append(Form("%s",sddstatus.Data()));
 cfcontname_casccuts.Append(Form("%s",fSuffix.Data()));
 // - CFContainer PID study Xi minus
 if(!fCFContCascadePIDXiMinus)  {
   const Int_t  lNbSteps      =  7 ;
   const Int_t  lNbVariables  =  3 ;
     //Array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[3] = {0};
   lNbBinsPerVar[0] = 100;
   lNbBinsPerVar[1] = 800;
   lNbBinsPerVar[2] = 22;
   fCFContCascadePIDXiMinus = new AliCFContainer(cfcontname_cascpidximinus,"Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
     //Setting the bin limits 
   fCFContCascadePIDXiMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
   fCFContCascadePIDXiMinus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
   fCFContCascadePIDXiMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
     //Setting the step title : one per PID case
   fCFContCascadePIDXiMinus->SetStepTitle(0, "No PID");
   fCFContCascadePIDXiMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
   fCFContCascadePIDXiMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
   fCFContCascadePIDXiMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
   fCFContCascadePIDXiMinus->SetStepTitle(4, "Comb. PID / Bachelor");
   fCFContCascadePIDXiMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
   fCFContCascadePIDXiMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");  
     //Setting the variable title, per axis
   fCFContCascadePIDXiMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
   fCFContCascadePIDXiMinus->SetVarTitle(1, "M( #Lambda , #pi^{-} ) (GeV/c^{2})");
   fCFContCascadePIDXiMinus->SetVarTitle(2, "Y_{#Xi}");
 }
 // - CFContainer PID study Xi plus
 if (!fCFContCascadePIDXiPlus) {
   const Int_t  lNbSteps      =  7 ;
   const Int_t  lNbVariables  =  3 ;
     //Array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[3] = {0};
   lNbBinsPerVar[0] = 100;
   lNbBinsPerVar[1] = 800;
   lNbBinsPerVar[2] = 22;
   fCFContCascadePIDXiPlus = new AliCFContainer(cfcontname_cascpidxiplus,"Pt_{cascade} Vs M_{#bar{#Xi}^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
     //Setting the bin limits 
   fCFContCascadePIDXiPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
   fCFContCascadePIDXiPlus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
   fCFContCascadePIDXiPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
     //Setting the step title : one per PID case
   fCFContCascadePIDXiPlus->SetStepTitle(0, "No PID");
   fCFContCascadePIDXiPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
   fCFContCascadePIDXiPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
   fCFContCascadePIDXiPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
   fCFContCascadePIDXiPlus->SetStepTitle(4, "Comb. PID / Bachelor");
   fCFContCascadePIDXiPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
   fCFContCascadePIDXiPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
     //Setting the variable title, per axis
   fCFContCascadePIDXiPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
   fCFContCascadePIDXiPlus->SetVarTitle(1, "M( #Lambda , #pi^{+} ) (GeV/c^{2})");
   fCFContCascadePIDXiPlus->SetVarTitle(2, "Y_{#Xi}");
 }
 // - CFContainer PID study Omega minus
 if(!fCFContCascadePIDOmegaMinus)  {
   const Int_t  lNbSteps      =  7 ;
   const Int_t  lNbVariables  =  3 ;
     //Array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[3] = {0};
   lNbBinsPerVar[0] = 100;
   lNbBinsPerVar[1] = 1000;
   lNbBinsPerVar[2] = 22;
   fCFContCascadePIDOmegaMinus = new AliCFContainer(cfcontname_cascpidomegaminus,"Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
     //Setting the bin limits 
   fCFContCascadePIDOmegaMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
   fCFContCascadePIDOmegaMinus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
   fCFContCascadePIDOmegaMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
     //Setting the step title : one per PID case
   fCFContCascadePIDOmegaMinus->SetStepTitle(0, "No PID");
   fCFContCascadePIDOmegaMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
   fCFContCascadePIDOmegaMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
   fCFContCascadePIDOmegaMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
   fCFContCascadePIDOmegaMinus->SetStepTitle(4, "Comb. PID / Bachelor");
   fCFContCascadePIDOmegaMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
   fCFContCascadePIDOmegaMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
     //Setting the variable title, per axis
   fCFContCascadePIDOmegaMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
   fCFContCascadePIDOmegaMinus->SetVarTitle(1, "M( #Lambda , K^{-} ) (GeV/c^{2})");
   fCFContCascadePIDOmegaMinus->SetVarTitle(2, "Y_{#Omega}");
 }
 // - CFContainer PID study Omega plus
 if(!fCFContCascadePIDOmegaPlus)  {
   const Int_t  lNbSteps      =  7 ;
   const Int_t  lNbVariables  =  3 ;
     //Array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[3] = {0};
   lNbBinsPerVar[0] = 100;
   lNbBinsPerVar[1] = 1000;
   lNbBinsPerVar[2] = 22; 
   fCFContCascadePIDOmegaPlus = new AliCFContainer(cfcontname_cascpidomegaplus,"Pt_{cascade} Vs M_{#bar{#Omega}^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
     //Setting the bin limits 
   fCFContCascadePIDOmegaPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
   fCFContCascadePIDOmegaPlus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
   fCFContCascadePIDOmegaPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity 
     //Setting the step title : one per PID case
   fCFContCascadePIDOmegaPlus->SetStepTitle(0, "No PID");
   fCFContCascadePIDOmegaPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
   fCFContCascadePIDOmegaPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
   fCFContCascadePIDOmegaPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
   fCFContCascadePIDOmegaPlus->SetStepTitle(4, "Comb. PID / Bachelor");
   fCFContCascadePIDOmegaPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
   fCFContCascadePIDOmegaPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
     //Setting the variable title, per axis
   fCFContCascadePIDOmegaPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
   fCFContCascadePIDOmegaPlus->SetVarTitle(1, "M( #Lambda , K^{+} ) (GeV/c^{2})");
   fCFContCascadePIDOmegaPlus->SetVarTitle(2, "Y_{#Omega}");  
 }
 // - CFContainer: towards the optimisation of topological selections
 if(! fCFContCascadeCuts) {
	// Container meant to store all the relevant distributions corresponding to the cut variables.
        // NB: overflow/underflow of variables on which we want to cut later should be 0!!! 
   const Int_t  lNbSteps      =  4 ;
   const Int_t  lNbVariables  =  20 ;
     //Array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[lNbVariables] = {0};
   lNbBinsPerVar[0]  = 25;     //DcaCascDaughters             :  [0.0,2.5]           -> Rec.Cut = 2.0;
   lNbBinsPerVar[1]  = 25;     //DcaBachToPrimVertex          :  [0.0,0.24,100.0]    -> Rec.Cut = 0.01; 
   lNbBinsPerVar[2]  = 61;     //CascCosineOfPointingAngle    :  [0.94,1.001]        -> Rec.Cut = 0.95;
   lNbBinsPerVar[3]  = 40;     //CascRadius                   :  [0.0,3.9,1000.0]    -> Rec.Cut = 0.2;
   lNbBinsPerVar[4]  = 30;     //InvMassLambdaAsCascDghter    :  [1.1,1.3]           -> Rec.Cut = 0.008;
   lNbBinsPerVar[5]  = 21;     //DcaV0Daughters               :  [0.0,2.0]           -> Rec.Cut = 1.5;
   lNbBinsPerVar[6]  = 201;    //V0CosineOfPointingAngleToXi  :  [0.89,1.0]          -> No Rec.Cut;
   lNbBinsPerVar[7]  = 40;     //V0Radius                     :  [0.0,3.9,1000.0]    -> Rec.Cut = 0.2;
   lNbBinsPerVar[8]  = 40;     //DcaV0ToPrimVertex            :  [0.0,0.39,110.0]    -> Rec.Cut = 0.01;  
   lNbBinsPerVar[9]  = 25;     //DcaPosToPrimVertex           :  [0.0,0.24,100.0]    -> Rec.Cut = 0.01;
   lNbBinsPerVar[10] = 25;     //DcaNegToPrimVertex           :  [0.0,0.24,100.0]    -> Rec.Cut = 0.01;
   lNbBinsPerVar[11] = 250;    //InvMassXi                    :   1-MeV/c2 bins
   lNbBinsPerVar[12] = 220;    //InvMassOmega                 :   1-MeV/c2 bins
   lNbBinsPerVar[13] = 100;    //XiTransvMom                  :  [0.0,10.0]
   lNbBinsPerVar[14] = 110;    //Y(Xi)                        :   0.02 in rapidity units
   lNbBinsPerVar[15] = 110;    //Y(Omega)                     :   0.02 in rapidity units
   lNbBinsPerVar[16] = 112;    //Proper lenght of cascade       
   lNbBinsPerVar[17] = 112;    //Proper lenght of V0
   lNbBinsPerVar[18] = 112;    //Distance V0-Xi in transverse plane
   lNbBinsPerVar[19] = 26;     //Proton cascade daughter momentum
   fCFContCascadeCuts = new AliCFContainer(cfcontname_casccuts,"Container for Cascade cuts", lNbSteps, lNbVariables, lNbBinsPerVar);
     //Setting the bin limits 
     //0 -  DcaXiDaughters
   //Double_t *lBinLim0  = new Double_t[ lNbBinsPerVar[0] + 1 ];
   //     for(Int_t i=0; i< lNbBinsPerVar[0]; i++) lBinLim0[i] = (Double_t)0.0 + (2.4 - 0.0)/(lNbBinsPerVar[0] - 1) * (Double_t)i;
   //     lBinLim0[ lNbBinsPerVar[0] ] = 3.0;
   //fCFContCascadeCuts -> SetBinLimits(0, lBinLim0);  
   //delete [] lBinLim0;
   fCFContCascadeCuts->SetBinLimits(0,0.0,2.5);
     //1 - DcaToPrimVertexXi
   Double_t *lBinLim1  = new Double_t[ lNbBinsPerVar[1] + 1 ];
        for(Int_t i=0; i<lNbBinsPerVar[1]; i++) lBinLim1[i] = (Double_t)0.0 + (0.24  - 0.0)/(lNbBinsPerVar[1] - 1) * (Double_t)i;
        lBinLim1[ lNbBinsPerVar[1] ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(1, lBinLim1);  
   delete [] lBinLim1;    
     //2 - CascCosineOfPointingAngle 
   fCFContCascadeCuts->SetBinLimits(2, 0.94, 1.001);
     //3 - CascRadius
   Double_t *lBinLim3  = new Double_t[ lNbBinsPerVar[3]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[3]; i++)   lBinLim3[i]  = (Double_t)0.0   + (3.9  - 0.0 )/(lNbBinsPerVar[3] - 1)  * (Double_t)i ;
        lBinLim3[ lNbBinsPerVar[3] ] = 1000.0;
   fCFContCascadeCuts -> SetBinLimits(3,  lBinLim3 );        
   delete [] lBinLim3;
     //4 - InvMassLambdaAsCascDghter
   fCFContCascadeCuts->SetBinLimits(4, 1.1, 1.13);
     //5 - DcaV0Daughters
   fCFContCascadeCuts -> SetBinLimits(5, 0., 2.1);
     //6 - V0CosineOfPointingAngle
   fCFContCascadeCuts -> SetBinLimits(6, 0.8, 1.001);
     //7 - V0Radius
   Double_t *lBinLim7 = new Double_t[ lNbBinsPerVar[7] + 1];
        for(Int_t i=0; i< lNbBinsPerVar[7];i++) lBinLim7[i] = (Double_t)0.0 + (3.9 - 0.0)/(lNbBinsPerVar[7] - 1) * (Double_t)i;
        lBinLim7[ lNbBinsPerVar[7] ] = 1000.0;
   fCFContCascadeCuts -> SetBinLimits(7, lBinLim7); 
   delete [] lBinLim7;
     //8 - DcaV0ToPrimVertex
   Double_t *lBinLim8  = new Double_t[ lNbBinsPerVar[8]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[8];i++)   lBinLim8[i]  = (Double_t)0.0   + (0.39  - 0.0 )/(lNbBinsPerVar[8]-1)  * (Double_t)i ;
        lBinLim8[ lNbBinsPerVar[8]  ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(8,  lBinLim8 );    
   delete [] lBinLim8;
     //9 - DcaPosToPrimVertex
   Double_t *lBinLim9  = new Double_t[ lNbBinsPerVar[9]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[9];i++)   lBinLim9[i]  = (Double_t)0.0   + (0.24  - 0.0 )/(lNbBinsPerVar[9]-1)  * (Double_t)i ;
        lBinLim9[ lNbBinsPerVar[9]  ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(9,  lBinLim9 );        
   delete [] lBinLim9;
     //10 - DcaNegToPrimVertex
   Double_t *lBinLim10  = new Double_t[ lNbBinsPerVar[10]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[10];i++)   lBinLim10[i]  = (Double_t)0.0   + (0.24  - 0.0 )/(lNbBinsPerVar[10]-1)  * (Double_t)i ;
        lBinLim10[ lNbBinsPerVar[10]  ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(10,  lBinLim10 );            // DcaPosToPrimVertexXi : 0.0 to 0.25 
   delete [] lBinLim10;
     //11 - InvMassXi
   fCFContCascadeCuts->SetBinLimits(11, 1.20, 1.45);
     //12 - InvMassOmega
   fCFContCascadeCuts->SetBinLimits(12, 1.57, 1.79);
     //13 - XiTransvMom
   fCFContCascadeCuts->SetBinLimits(13, 0.0, 10.0); 
     //14 - Y(Xi)
   fCFContCascadeCuts->SetBinLimits(14, -1.1, 1.1);
     //15 - Y(Omega)
   fCFContCascadeCuts->SetBinLimits(15, -1.1, 1.1);
     //16 - Proper time of cascade
   Double_t *lBinLim16  = new Double_t[ lNbBinsPerVar[16]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[16];i++) lBinLim16[i] = (Double_t) -1. + (110. + 1.0 ) / (lNbBinsPerVar[16] - 1) * (Double_t) i;
   lBinLim16[ lNbBinsPerVar[16] ] = 2000.0;
   fCFContCascadeCuts->SetBinLimits(16, lBinLim16);
     //17 - Proper time of V0
   fCFContCascadeCuts->SetBinLimits(17, lBinLim16);
     //18 - Distance V0-Xi in transverse plane
   fCFContCascadeCuts->SetBinLimits(18, lBinLim16);
     //19 - (Anti-)Proton cascade daughter momentum
   Double_t *lBinLim19  = new Double_t[ lNbBinsPerVar[19]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[19];i++)   lBinLim19[i] = (Double_t)0.0 + (0.5-0.0)/(lNbBinsPerVar[19]-1) * (Double_t)i;
        lBinLim19[ lNbBinsPerVar[19]  ] = 50.0;
   fCFContCascadeCuts->SetBinLimits(19, lBinLim19);
   delete [] lBinLim19;
     // Setting the number of steps : one for each cascade species (Xi-, Xi+ and Omega-, Omega+)
   fCFContCascadeCuts->SetStepTitle(0, "#Xi^{-} candidates");
   fCFContCascadeCuts->SetStepTitle(1, "#bar{#Xi}^{+} candidates");
   fCFContCascadeCuts->SetStepTitle(2, "#Omega^{-} candidates");
   fCFContCascadeCuts->SetStepTitle(3, "#bar{#Omega}^{+} candidates");
     // Setting the variable title, per axis
   fCFContCascadeCuts->SetVarTitle(0,  "Dca(cascade daughters) (cm)");
   fCFContCascadeCuts->SetVarTitle(1,  "ImpactParamToPV(bachelor) (cm)");
   fCFContCascadeCuts->SetVarTitle(2,  "cos(cascade PA)");
   fCFContCascadeCuts->SetVarTitle(3,  "R_{2d}(cascade decay) (cm)");
   fCFContCascadeCuts->SetVarTitle(4,  "M_{#Lambda}(as casc dghter) (GeV/c^{2})");
   fCFContCascadeCuts->SetVarTitle(5,  "Dca(V0 daughters) in Xi (cm)");
   if      (fCollidingSystem == 0) fCFContCascadeCuts->SetVarTitle(6,  "cos(V0 PA) to cascade vtx");
   else if (fCollidingSystem == 1) fCFContCascadeCuts->SetVarTitle(6,  "cos(V0 PA) to primary vtx");
   fCFContCascadeCuts->SetVarTitle(7,  "R_{2d}(V0 decay) (cm)");
   fCFContCascadeCuts->SetVarTitle(8,  "ImpactParamToPV(V0) (cm)");
   fCFContCascadeCuts->SetVarTitle(9,  "ImpactParamToPV(Pos) (cm)");
   fCFContCascadeCuts->SetVarTitle(10, "ImpactParamToPV(Neg) (cm)");
   fCFContCascadeCuts->SetVarTitle(11, "Inv. Mass(Xi) (GeV/c^{2})");
   fCFContCascadeCuts->SetVarTitle(12, "Inv. Mass(Omega) (GeV/c^{2})");
   fCFContCascadeCuts->SetVarTitle(13, "pt(cascade) (GeV/c)");
   fCFContCascadeCuts->SetVarTitle(14, "Y(Xi)");
   fCFContCascadeCuts->SetVarTitle(15, "Y(Omega)");
   fCFContCascadeCuts->SetVarTitle(16, "mL/p (cascade) (cm)");
   fCFContCascadeCuts->SetVarTitle(17, "mL/p (V0) (cm)");
   fCFContCascadeCuts->SetVarTitle(18, "Distance V0-Cascade in transverse plane (cm)");
   fCFContCascadeCuts->SetVarTitle(19, "(Anti-)Proton casc. daught. momemtum (GeV/c)");
 }

 PostData(1, fListHistCascade);
 PostData(2, fCFContCascadePIDXiMinus);
 PostData(3, fCFContCascadePIDXiPlus);
 PostData(4, fCFContCascadePIDOmegaMinus);
 PostData(5, fCFContCascadePIDOmegaPlus);
 PostData(6, fCFContCascadeCuts);
} // end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskCheckCascadepp::UserExec(Option_t *) {


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Main loop (called for each event)
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  //----------------
  //Define variables 
  AliESDEvent *lESDevent = 0x0;
  AliAODEvent *lAODevent = 0x0;

  //----------------------
  // Check the PIDresponse
  //---------------------
  if(!fPIDResponse) {
       AliError("Cannot get pid response");
       return;
  }
  if(! fESDtrackCuts ){
       fESDtrackCuts = new AliESDtrackCuts();
  }
  if(! fUtils ){
       fUtils = new AliAnalysisUtils();
  }

  ///////////////////
  // EVENT SELECTIONS
  ///////////////////
  // In order:
  // 0) SDD status 
  // 1) Pre-Trigger selections
  //     1.1) Incomplete DAQ events (introduced for Run2 2015 data)
  //     1.2) Background rejection based on SPD cluster vs tracklet correlation
  //     1.3) Pileup
  // 2) Trigger selection (Physics selection)
  // -) Cascade and V0 re-vertexer
  // 3) Well-established PV
  //     3.1) not only TPC vertex
  //     3.2) requirement on the resolution and dispersion
  //     3.2) distance between the two vertices
  //     3.4) |Zpv| < 10 cm
  // - Define useful variables

  // - Define useful variables
  Int_t ncascades = 0;
  Int_t nTrackMultiplicity = 0;

  //============================
  // Plots before any selections
  //============================
  if (fAnalysisType == "ESD") {
      lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
      if (!lESDevent) {
         AliWarning("ERROR: lESDevent not available \n");
         return;
      }
      ncascades          = lESDevent->GetNumberOfCascades();
      nTrackMultiplicity = fESDtrackCuts->GetReferenceMultiplicity(lESDevent,AliESDtrackCuts::kTrackletsITSTPC,0.5);
  } else if (fAnalysisType == "AOD") {
      lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
      if (!lAODevent) {
          AliWarning("ERROR: lAODevent not available \n");
          return;
      }
      ncascades          = lAODevent->GetNumberOfCascades();
      nTrackMultiplicity = -100;  //FIXME: I can't find the equivalent method for the AOD  
  } else {
      Printf("Analysis type (ESD or AOD) not specified \n");
      return;
  }
  fHistCascadeMultiplicityBeforeAnySel->Fill(ncascades);

  //========================
  // 0) SDD status selection
  //========================
  if (fApplyEvSelSDDstatus && fCollidingSystem == 0) {
       if (fAnalysisType == "ESD") {
           TString trcl = lESDevent->GetFiredTriggerClasses();
           if (fwithSDD && !(trcl.Contains("ALLNOTRD"))) {
               AliWarning("We are selecting events with SDD turn ON. This event has the SDD turn OFF. =>  RETURN!! (Exclude it)...");
               PostData(1, fListHistCascade);
               PostData(2, fCFContCascadePIDXiMinus);
               PostData(3, fCFContCascadePIDXiPlus);
               PostData(4, fCFContCascadePIDOmegaMinus);
               PostData(5, fCFContCascadePIDOmegaPlus);
               PostData(6, fCFContCascadeCuts);
               return;
           } else if (!fwithSDD && (trcl.Contains("ALLNOTRD"))) {
               AliWarning("We are selecting events with SDD turn OFF. This event has the SDD turn ON. =>  RETURN!! (Exclude it)...");
               PostData(1, fListHistCascade);
               PostData(2, fCFContCascadePIDXiMinus);
               PostData(3, fCFContCascadePIDXiPlus);
               PostData(4, fCFContCascadePIDOmegaMinus);
               PostData(5, fCFContCascadePIDOmegaPlus);
               PostData(6, fCFContCascadeCuts);
               return;
           }
       } else if (fAnalysisType == "AOD") {
           TString trcl = lAODevent->GetFiredTriggerClasses();
           if (fwithSDD && !(trcl.Contains("ALLNOTRD"))) {
               AliWarning("We are selecting events with SDD turn ON. This event has the SDD turn OFF. =>  RETURN!! (Exclude it)...");
               PostData(1, fListHistCascade);
               PostData(2, fCFContCascadePIDXiMinus);
               PostData(3, fCFContCascadePIDXiPlus);
               PostData(4, fCFContCascadePIDOmegaMinus);
               PostData(5, fCFContCascadePIDOmegaPlus);
               PostData(6, fCFContCascadeCuts);
               return;
           } else if (!fwithSDD && (trcl.Contains("ALLNOTRD"))) {
               AliWarning("We are selecting events with SDD turn OFF. This event has the SDD turn ON. =>  RETURN!! (Exclude it)...");
               PostData(1, fListHistCascade);
               PostData(2, fCFContCascadePIDXiMinus);
               PostData(3, fCFContCascadePIDXiPlus);
               PostData(4, fCFContCascadePIDOmegaMinus);
               PostData(5, fCFContCascadePIDOmegaPlus);
               PostData(6, fCFContCascadeCuts);
               return;
           }
       }
  }
  fHistCascadeMultiplicityAfterSDDstatusSel->Fill(ncascades);

  //======================================
  // 1.1) Removal of incomplete DAQ events
  //======================================
  Bool_t IncompleteDAQ = kFALSE;
  if (fAnalysisType == "ESD") IncompleteDAQ = lESDevent->IsIncompleteDAQ();
  if (fApplyEvSelDAQincomplete && IncompleteDAQ){
        AliWarning("This is a DAQ incomplete event... return !");
        PostData(1, fListHistCascade);
        PostData(2, fCFContCascadePIDXiMinus);
        PostData(3, fCFContCascadePIDXiPlus);
        PostData(4, fCFContCascadePIDOmegaMinus);
        PostData(5, fCFContCascadePIDOmegaPlus);
        PostData(6, fCFContCascadeCuts);
        return;
  }
  fHistCascadeMultiplicityAfterDAQincompleteEvRej->Fill(ncascades);

  //=======================================================================
  // 1.2) Background rejection based on SPD cluster vs tracklet correlation
  //=======================================================================
  if(fApplyEvSelSPDclustervstracklet){
      if (fAnalysisType == "ESD" && fUtils->IsSPDClusterVsTrackletBG(lESDevent)) {
              AliWarning("Pb / Is background based on SPD cluster vs tracklet correlation... return!");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDXiMinus);
              PostData(3, fCFContCascadePIDXiPlus);
              PostData(4, fCFContCascadePIDOmegaMinus);
              PostData(5, fCFContCascadePIDOmegaPlus);
              PostData(6, fCFContCascadeCuts);
              return;
      }
  }
  fHistCascadeMultiplicityAfterSPDclustervstrackletSel->Fill(ncascades);

  //======================
  // 1.3) Pileup selection
  //======================
  if (fApplyEvSelPileup && fCollidingSystem == 0) {
      if (fAnalysisType == "ESD") {
           if(lESDevent->IsPileupFromSPD(fSPDPileUpminContributors)){
              AliWarning("Pb / Pile-up event ... return!");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDXiMinus);
              PostData(3, fCFContCascadePIDXiPlus);
              PostData(4, fCFContCascadePIDOmegaMinus);
              PostData(5, fCFContCascadePIDOmegaPlus);
              PostData(6, fCFContCascadeCuts);
              return;
           }
      } else if (fAnalysisType == "AOD") {
           if(lAODevent->IsPileupFromSPD(fSPDPileUpminContributors)){
              AliWarning("Pb / Pile-up event ... return!");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDXiMinus);
              PostData(3, fCFContCascadePIDXiPlus);
              PostData(4, fCFContCascadePIDOmegaMinus);
              PostData(5, fCFContCascadePIDOmegaPlus);
              PostData(6, fCFContCascadeCuts);
              return;
           }
      }
  }
  fHistCascadeMultiplicityAfterPileupRej->Fill(ncascades);

  //========================
  // Multiplicity estimators   //FIXME: not implemented at the moment
  //========================
  // For the pPb collisions
  //AliCentrality* centrality = 0;
  //if      (fAnalysisType == "ESD" && (fCollidingSystem == 1)) centrality = lESDevent->GetCentrality();
  //else if (fAnalysisType == "AOD" && (fCollidingSystem == 1)) centrality = lAODevent->GetCentrality();
  //Float_t lcentrality = 0.;
  //if (fCollidingSystem == 1) {
  //     if (fkUseCleaning) lcentrality = centrality->GetCentralityPercentile(fCentrEstimator.Data());
  //     else {
  //         lcentrality = centrality->GetCentralityPercentileUnchecked(fCentrEstimator.Data());
  //         if (centrality->GetQuality()>1) {
  //             PostData(1, fCFContCascadeCuts);
  //            return;
  //         }
  //     }
  //     //if (lcentrality<fCentrLowLim||lcentrality>=fCentrUpLim) { 
  //     //    PostData(1, fCFContCascadeCuts);
  //     //    return;
  //     //}
  //} else if (fCollidingSystem == 0) lcentrality = 0.;  

  //=====================
  // 2) Physics selection 
  //=====================
  if (fApplyEvSelPhysicsSel) {
      UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
      Bool_t isSelected = 0;
      isSelected = (maskIsSelected & fkTriggerClass) == fkTriggerClass;
      if (!isSelected){
          AliWarning("We are selecting the events that past the Physics Selection. This event does not pass the Physics Selection. =>  RETURN!! (Exclude it)..."); 
          PostData(1, fListHistCascade);
          PostData(2, fCFContCascadePIDXiMinus);
          PostData(3, fCFContCascadePIDXiPlus);
          PostData(4, fCFContCascadePIDOmegaMinus);
          PostData(5, fCFContCascadePIDOmegaPlus);
          PostData(6, fCFContCascadeCuts);
          return;
      }
  }
  fHistCascadeMultiplicityAfterPhysicsSel->Fill(ncascades);

  //=========================================
  // - Re-run cascade vertexer (only for ESD)
  //=========================================
  if (fAnalysisType == "ESD" && fRerunV0CascVertexers) {
       lESDevent->ResetCascades();
       lESDevent->ResetV0s();
       AliV0vertexer *lV0vtxer = new AliV0vertexer();
       AliCascadeVertexer *lCascVtxer = new AliCascadeVertexer();
       //lV0vtxer->GetCuts(fV0Sels);
       //lCascVtxer->GetCuts(fCascSels);
       lV0vtxer->SetCuts(fV0Sels);      // NB don't use SetDefaultCuts!! because it acts on static variables 
       lCascVtxer->SetCuts(fCascSels);
       lV0vtxer->Tracks2V0vertices(lESDevent);
       lCascVtxer->V0sTracks2CascadeVertices(lESDevent);
       //delete lV0vtxer;
       //delete lCascVtxer;
  }
  fHistCascadeMultiplicityAfterRevertexing->Fill(ncascades);

  //====================================
  // 3) Primary Vertex quality selection
  //====================================
  const AliESDVertex *lESDPrimaryTrackingVtx = 0x0;
  const AliESDVertex *lESDPrimarySPDVtx      = 0x0;
  const AliAODVertex *lAODPrimaryTrackingVtx = 0x0;
  const AliAODVertex *lAODPrimarySPDVtx      = 0x0;
  if (fAnalysisType == "ESD") {
       lESDPrimaryTrackingVtx = lESDevent->GetPrimaryVertexTracks();
       lESDPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
  } else if (fAnalysisType == "AOD") {
       lAODPrimarySPDVtx = lAODevent->GetPrimaryVertexSPD();
       lAODPrimaryTrackingVtx = lAODevent->GetPrimaryVertex(); 
  }

  // 3.1) reject events if both are explicitly requested and none is available
  if (fCollidingSystem == 0) {
       if (fApplyEvSelNoTPConlyPrimVtx) {
            if (!(lESDPrimarySPDVtx->GetStatus() && lESDPrimaryTrackingVtx->GetStatus()) && fAnalysisType == "ESD"){
                  AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDXiMinus);
                  PostData(3, fCFContCascadePIDXiPlus);
                  PostData(4, fCFContCascadePIDOmegaMinus);
                  PostData(5, fCFContCascadePIDOmegaPlus);
                  PostData(6, fCFContCascadeCuts);
                  return;
            }
            if (!(lAODPrimarySPDVtx && lAODPrimaryTrackingVtx) && fAnalysisType == "AOD") {
                  AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDXiMinus);
                  PostData(3, fCFContCascadePIDXiPlus);
                  PostData(4, fCFContCascadePIDOmegaMinus);
                  PostData(5, fCFContCascadePIDOmegaPlus);
                  PostData(6, fCFContCascadeCuts);
                  return;
             }
        }
  }
  fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel->Fill(ncascades);

  // 3.2) check the spd vertex resolution and reject if not satisfied //FIXME: only for ESD
  if (fCollidingSystem == 0) {
        if (!lESDPrimaryTrackingVtx->GetStatus()) {
             if (!lESDPrimarySPDVtx->GetStatus()) {
                  AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDXiMinus);
                  PostData(3, fCFContCascadePIDXiPlus);
                  PostData(4, fCFContCascadePIDOmegaMinus);
                  PostData(5, fCFContCascadePIDOmegaPlus);
                  PostData(6, fCFContCascadeCuts);
                  return;
             }
             if (fApplyEvSelSPDvtxres && lESDPrimarySPDVtx->GetStatus() && lESDPrimarySPDVtx->IsFromVertexerZ() && !(lESDPrimarySPDVtx->GetDispersion()<0.04 && lESDPrimarySPDVtx->GetZRes()<0.25)) {
                  AliWarning("Pb / The SPD prim. vertex has a Z resolution > 0.25 and dispersion > 0.04 ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDXiMinus);
                  PostData(3, fCFContCascadePIDXiPlus);
                  PostData(4, fCFContCascadePIDOmegaMinus);
                  PostData(5, fCFContCascadePIDOmegaPlus);
                  PostData(6, fCFContCascadeCuts);
                  return;
             }
        } else {
             if (fApplyEvSelSPDvtxres && lESDPrimarySPDVtx->GetStatus() && lESDPrimarySPDVtx->IsFromVertexerZ() && !(lESDPrimarySPDVtx->GetDispersion()<0.04 && lESDPrimarySPDVtx->GetZRes()<0.25)) {
                  AliWarning("Pb / The SPD prim. vertex has a Z resolution > 0.25 and dispersion > 0.04 ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDXiMinus);
                  PostData(3, fCFContCascadePIDXiPlus);
                  PostData(4, fCFContCascadePIDOmegaMinus);
                  PostData(5, fCFContCascadePIDOmegaPlus);
                  PostData(6, fCFContCascadeCuts);
                  return;
             }
        }
  }
  fHistCascadeMultiplicityAfterSPDresolution->Fill(ncascades);

  // 3.3) check the proximity between the spd vertex and trak vertex, and reject if not satisfied 
  if (fCollidingSystem == 0) {
        if (lESDPrimaryTrackingVtx->GetStatus() && lESDPrimarySPDVtx->GetStatus()) {
             if (fApplyEvSelVtxProximity && (TMath::Abs(lESDPrimarySPDVtx->GetZ() - lESDPrimaryTrackingVtx->GetZ()) > 0.5)) {
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDXiMinus);
                  PostData(3, fCFContCascadePIDXiPlus);
                  PostData(4, fCFContCascadePIDOmegaMinus);
                  PostData(5, fCFContCascadePIDOmegaPlus);
                  PostData(6, fCFContCascadeCuts);
                  return;
             }
        }
  } else if (fCollidingSystem == 1) { //FIXME: the "activation" variable ares not used
          if (fAnalysisType == "ESD") {
              Bool_t fHasVertex = kFALSE;
              const AliESDVertex *vertex = lESDevent->GetPrimaryVertexTracks();
              if (vertex->GetNContributors() < 1) {
                  vertex = lESDevent->GetPrimaryVertexSPD();
                  if (vertex->GetNContributors() < 1) fHasVertex = kFALSE;
                  else                                fHasVertex = kTRUE;
                  Double_t cov[6]={0};
                  vertex->GetCovarianceMatrix(cov);
                  Double_t zRes = TMath::Sqrt(cov[5]);
                  if (vertex->IsFromVertexerZ() && (zRes>0.25)) fHasVertex = kFALSE;
              } else fHasVertex = kTRUE;
              if (fHasVertex == kFALSE) {
                  AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDXiMinus);
                  PostData(3, fCFContCascadePIDXiPlus);
                  PostData(4, fCFContCascadePIDOmegaMinus);
                  PostData(5, fCFContCascadePIDOmegaPlus);
                  PostData(6, fCFContCascadeCuts);
                  return;
              }
              if (fUtils->IsFirstEventInChunk(lESDevent)) { //Is First event in chunk rejection: Still present!
                  AliWarning("Pb / This is the first event in the chunk! ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDXiMinus);
                  PostData(3, fCFContCascadePIDXiPlus);
                  PostData(4, fCFContCascadePIDOmegaMinus);
                  PostData(5, fCFContCascadePIDOmegaPlus);
                  PostData(6, fCFContCascadeCuts);
                  return;
              }
          } else if (fAnalysisType == "AOD") {
              Bool_t fHasVertex = kFALSE;
              const AliAODVertex *vertex = lAODevent->GetPrimaryVertex();
              if (vertex->GetNContributors() < 1) {
                  vertex = lAODevent->GetPrimaryVertexSPD();
                  if (vertex->GetNContributors() < 1) fHasVertex = kFALSE;
                  else                                fHasVertex = kTRUE;
                  Double_t cov[6]={0};
                  vertex->GetCovarianceMatrix(cov);
                  Double_t zRes = TMath::Sqrt(cov[5]);
                  if (vertex->IsFromVertexerZ() && (zRes>0.25)) fHasVertex = kFALSE;
              } else fHasVertex = kTRUE;
              if (fHasVertex == kFALSE) {
                  AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDXiMinus);
                  PostData(3, fCFContCascadePIDXiPlus);
                  PostData(4, fCFContCascadePIDOmegaMinus);
                  PostData(5, fCFContCascadePIDOmegaPlus);
                  PostData(6, fCFContCascadeCuts);
                  return;
              }
              if (fHasVertex == kFALSE) { //Is First event in chunk rejection: Still present!  //FIXME
                  AliWarning("Pb / This is the first event in the chunk! ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDXiMinus);
                  PostData(3, fCFContCascadePIDXiPlus);
                  PostData(4, fCFContCascadePIDOmegaMinus);
                  PostData(5, fCFContCascadePIDOmegaPlus);
                  PostData(6, fCFContCascadeCuts);
                  return;
              }  
          }
      }
  fHistCascadeMultiplicityAfterVerticesProximity->Fill(ncascades);
 
  //=========================================================
  // 3.4) Vertex Z position selection (+ magnetic field info)
  //=========================================================
  Double_t lBestPrimaryVtxPos[3]  = {-100.0, -100.0, -100.0};
  Double_t tPrimaryVtxPosition[3] = {-100.0, -100.0, -100.0};
  if (fAnalysisType == "ESD") {
      const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex(); 
      if (!lPrimaryBestESDVtx){
          AliWarning("No prim. vertex in ESD... return!");
          PostData(1, fListHistCascade);
          PostData(2, fCFContCascadePIDXiMinus);
          PostData(3, fCFContCascadePIDXiPlus);
          PostData(4, fCFContCascadePIDOmegaMinus);
          PostData(5, fCFContCascadePIDOmegaPlus);
          PostData(6, fCFContCascadeCuts);
          return;
      }
      lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
      // - Fill the vertex plots before any event selection on vertex position
      const AliVVertex *primaryVtx = lESDevent->GetPrimaryVertex();
      tPrimaryVtxPosition[0] = primaryVtx->GetX();
      tPrimaryVtxPosition[1] = primaryVtx->GetY();
      tPrimaryVtxPosition[2] = primaryVtx->GetZ();
  } else if (fAnalysisType == "AOD") {
      const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();
      if (!lPrimaryBestAODVtx){
          AliWarning("No prim. vertex in AOD... return!");
          PostData(1, fListHistCascade);
          PostData(2, fCFContCascadePIDXiMinus);
          PostData(3, fCFContCascadePIDXiPlus);
          PostData(4, fCFContCascadePIDOmegaMinus);
          PostData(5, fCFContCascadePIDOmegaPlus);
          PostData(6, fCFContCascadeCuts);
          return;
      }
      lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );
      // - Fill the vertex plots before any event selection on vertex position
      const AliVVertex *primaryVtx = lAODevent->GetPrimaryVertex();
      tPrimaryVtxPosition[0] = primaryVtx->GetX();
      tPrimaryVtxPosition[1] = primaryVtx->GetY();
      tPrimaryVtxPosition[2] = primaryVtx->GetZ();
  }
  fHistPVx->Fill( tPrimaryVtxPosition[0] );
  fHistPVy->Fill( tPrimaryVtxPosition[1] );
  fHistPVz->Fill( tPrimaryVtxPosition[2] );
  // - Selection on the primary vertex Z position  
  if (fApplyEvSelZprimVtxPos) {
      if (fAnalysisType == "ESD") {
          if (TMath::Abs(lBestPrimaryVtxPos[2]) < fVtxRangeMin || TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRangeMax) {
               AliWarning("Pb / | Z position of Best Prim Vtx | out of range [Rmin, Rmax]... return !");
               PostData(1, fListHistCascade);
               PostData(2, fCFContCascadePIDXiMinus);
               PostData(3, fCFContCascadePIDXiPlus);
               PostData(4, fCFContCascadePIDOmegaMinus);
               PostData(5, fCFContCascadePIDOmegaPlus);
               PostData(6, fCFContCascadeCuts);
               return;
          }
      } else if (fAnalysisType == "AOD") {
          if (TMath::Abs(lBestPrimaryVtxPos[2]) < fVtxRangeMin || TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRangeMax) {
              AliWarning("Pb / | Z position of Best Prim Vtx | out of range [Rmin, Rmax]... return !");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDXiMinus);
              PostData(3, fCFContCascadePIDXiPlus);
              PostData(4, fCFContCascadePIDOmegaMinus);
              PostData(5, fCFContCascadePIDOmegaPlus);
              PostData(6, fCFContCascadeCuts);
              return;
          }
      }
  }
  fHistCascadeMultiplicityAfterZprimVtxPosSel->Fill(ncascades);
  // - Vertex position plots: after any event selections
  tPrimaryVtxPosition[0] = -100.0;
  tPrimaryVtxPosition[1] = -100.0;
  tPrimaryVtxPosition[2] = -100.0;
  if (fAnalysisType == "ESD") {
       const AliVVertex *primaryVtx = lESDevent->GetPrimaryVertex();
       tPrimaryVtxPosition[0] = primaryVtx->GetX();
       tPrimaryVtxPosition[1] = primaryVtx->GetY();
       tPrimaryVtxPosition[2] = primaryVtx->GetZ();
  } else if (fAnalysisType == "AOD") {
       const AliVVertex *primaryVtx = lAODevent->GetPrimaryVertex();
       tPrimaryVtxPosition[0] = primaryVtx->GetX();
       tPrimaryVtxPosition[1] = primaryVtx->GetY();
       tPrimaryVtxPosition[2] = primaryVtx->GetZ();
  }
  fHistPVxAnalysis->Fill( tPrimaryVtxPosition[0] );
  fHistPVyAnalysis->Fill( tPrimaryVtxPosition[1] );
  fHistPVzAnalysis->Fill( tPrimaryVtxPosition[2] );
  
  // - Magnetic field
  Double_t lMagneticField = -10.;
  if      (fAnalysisType == "ESD") lMagneticField = lESDevent->GetMagneticField();
  else if (fAnalysisType == "AOD") lMagneticField = lAODevent->GetMagneticField();
  //if(TMath::Abs(lMagneticField ) < 10e-6) continue;


  //////////////////////////////
  // CASCADE RECONSTRUCTION PART
  //////////////////////////////
  
  //%%%%%%%%%%%%%
  // Cascade loop
  Int_t numcascades = 0;
  if      (fAnalysisType == "ESD") numcascades = lESDevent->GetNumberOfCascades();
  else if (fAnalysisType == "AOD") numcascades = lAODevent->GetNumberOfCascades();

  for (Int_t iXi = 0; iXi < numcascades; iXi++) {// This is the begining of the Cascade loop (ESD or AOD)
	   
      // -----------------------------------------------------------------------
      // - Initialisation of the local variables that will be needed for ESD/AOD

      // - 0th part of initialisation : around primary vertex ...
      //Double_t lBestPrimaryVtxRadius3D     = -500.0;
      // - 1st part of initialisation : variables needed to store AliESDCascade data members
      Double_t lEffMassXi                  = 0.;
      Double_t lDcaXiDaughters             = -1.;
      Double_t lXiCosineOfPointingAngle    = -1.;
      Double_t lPosXi[3]                   = { -1000.0, -1000.0, -1000.0 };
      Double_t lXiRadius                   = -1000. ;
      // - 2nd part of initialisation : Nbr of clusters within TPC for the 3 daughter cascade tracks
      Double_t lPosTPCClusters             = -1; // For ESD only ...
      Double_t lNegTPCClusters             = -1; // For ESD only ...
      Double_t lBachTPCClusters            = -1; // For ESD only ...
      Double_t lPosTPCFindClusters         = -1; // For ESD only ...
      Double_t lNegTPCFindClusters         = -1; // For ESD only ...
      Double_t lBachTPCFindClusters        = -1; // For ESD only ...
      Double_t lInnerWallMomCascDghters[3] = {-100., -100., -100.};
      Double_t lTPCSignalCascDghters   [3] = {-100., -100., -100.};
      // - 3rd part of initialisation : about V0 part in cascades
      Double_t lInvMassLambdaAsCascDghter  = 0.;
      Double_t lDcaV0DaughtersXi           = -1.;
      Double_t lDcaBachToPrimVertexXi      = -1.;
      Double_t lDcaV0ToPrimVertexXi        = -1.;
      Double_t lDcaPosToPrimVertexXi       = -1.;
      Double_t lDcaNegToPrimVertexXi       = -1.;
      Double_t lV0CosineOfPointingAngleXi  = -1. ;
      Double_t lPosV0Xi[3]                 = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
      Double_t lV0RadiusXi                 = -1000.0;
      Double_t lV0quality                  = 0.;
      // - 4th part of initialisation : Effective masses
      Double_t lInvMassXiMinus             = 0.;
      Double_t lInvMassXiPlus              = 0.;
      Double_t lInvMassOmegaMinus          = 0.;
      Double_t lInvMassOmegaPlus           = 0.;
      // - 5th part of initialisation : PID treatment
      Bool_t   lIsPosInXiProton            = kFALSE;
      Bool_t   lIsPosInXiPion              = kFALSE;
      Bool_t   lIsPosInOmegaProton         = kFALSE;
      Bool_t   lIsPosInOmegaPion           = kFALSE;
      Bool_t   lIsNegInXiProton            = kFALSE;
      Bool_t   lIsNegInXiPion              = kFALSE;
      Bool_t   lIsNegInOmegaProton         = kFALSE;
      Bool_t   lIsNegInOmegaPion           = kFALSE;
      Bool_t   lIsBachelorKaon             = kFALSE;
      Bool_t   lIsBachelorPion             = kFALSE; 
      Bool_t   lIsBachelorKaonForTPC       = kFALSE; 
      Bool_t   lIsBachelorPionForTPC       = kFALSE; 
      Bool_t   lIsNegPionForTPC            = kFALSE; 
      Bool_t   lIsPosPionForTPC            = kFALSE; 
      Bool_t   lIsNegProtonForTPC          = kFALSE; 
      Bool_t   lIsPosProtonForTPC          = kFALSE; 
      // - 6th part of initialisation : extra info for QA
      Double_t lXiMomX         = 0.;
      Double_t lXiMomY         = 0.;
      Double_t lXiMomZ         = 0.;
      Double_t lXiTransvMom    = 0.;
      Double_t lXiTotMom       = 0.;
      Double_t lV0PMomX        = 0.;
      Double_t lV0PMomY        = 0.;
      Double_t lV0PMomZ        = 0.;
      Double_t lV0NMomX        = 0.;
      Double_t lV0NMomY        = 0.;
      Double_t lV0NMomZ        = 0.;
      Double_t lV0TotMom       = 0.;
      Double_t lBachMomX       = 0.;
      Double_t lBachMomY       = 0.;
      Double_t lBachMomZ       = 0.;
      Double_t lBachTransvMom  = 0.;
      Double_t lBachTotMom     = 0.;
      Double_t lpTrackTransvMom  = 0.;
      Double_t lnTrackTransvMom  = 0.;
      Short_t  lChargeXi       = -2;
      Double_t lV0toXiCosineOfPointingAngle = 0.;
      Double_t lRapXi   = -20.0, lRapOmega = -20.0, lEta = -20.0, lTheta = 360., lPhi = 720.;
      Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
      Double_t etaPos = 0., etaNeg = 0., etaBach = 0.;
  
      if (fAnalysisType == "ESD") { 
  
          // -------------------------------------
          // - Load the cascades from the handler
          AliESDcascade *xi = lESDevent->GetCascade(iXi);
          if (!xi) continue;

          //----------------------------------------------------------------------------        
          // - Assigning the necessary variables for specific AliESDcascade data members	
          lV0quality = 0.;
          xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis: cascade = Xi-decay
          lEffMassXi  	           = xi->GetEffMassXi();
          lDcaXiDaughters	   = xi->GetDcaXiDaughters();
          lXiCosineOfPointingAngle = xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
	                                                               //Take care : the best available vertex should be used (like in AliCascadeVertexer)
          xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] ); 
          lXiRadius        	   = TMath::Sqrt( lPosXi[0]*lPosXi[0] + lPosXi[1]*lPosXi[1] );
		
          //-------------------------------------------------------------------------------------------------------------------------------
	  // - Around the tracks: Bach + V0 (ESD). Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
          UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
          UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
          UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
                                      //Care track label can be negative in MC production (linked with the track quality)
                                      //However = normally, not the case for track index ...
          // - Rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
          if (lBachIdx == lIdxNegXi) { AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue; }
          if (lBachIdx == lIdxPosXi) { AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue; }
          // - Get the track for the daughters
          AliESDtrack *pTrackXi	   = lESDevent->GetTrack( lIdxPosXi );
          AliESDtrack *nTrackXi	   = lESDevent->GetTrack( lIdxNegXi );
          AliESDtrack *bachTrackXi = lESDevent->GetTrack( lBachIdx );
          if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
               AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
               continue;
          }
          // - Get the number of TPC clusters for the daughters
          lPosTPCClusters   = pTrackXi->GetTPCClusterInfo(2,1);             // Old for #TPCclusters: ->GetTPCNcls();
          lNegTPCClusters   = nTrackXi->GetTPCClusterInfo(2,1);             // Old for #TPCclusters: ->GetTPCNcls();
          lBachTPCClusters  = bachTrackXi->GetTPCClusterInfo(2,1);          // Old for #TPCclusters: ->GetTPCNcls();
          // - Get the number TPC findable clusters for daughters
          lPosTPCFindClusters   = pTrackXi->GetTPCNclsF();                  // New
          lNegTPCFindClusters   = nTrackXi->GetTPCNclsF();                  // New
          lBachTPCFindClusters  = bachTrackXi->GetTPCNclsF();               // New

          //-------------------------------------
          // - Rejection of a poor quality tracks
          if (fTrackQualityCutTPCrefit) {
                // - Poor quality related to TPCrefit
                ULong_t pStatus    = pTrackXi->GetStatus();
                ULong_t nStatus    = nTrackXi->GetStatus();
                ULong_t bachStatus = bachTrackXi->GetStatus();
                if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
                if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
                if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
          }
          if (fTrackQualityCutnTPCcls) {
                // - Poor quality related to TPC clusters
                if (lPosTPCClusters  < fMinnTPCcls) { AliWarning(Form("Pb / V0 Pos. track has less than %i TPC clusters ... continue!",fMinnTPCcls)); continue; }
                if (lNegTPCClusters  < fMinnTPCcls) { AliWarning(Form("Pb / V0 Neg. track has less than %i TPC clusters ... continue!",fMinnTPCcls)); continue; }
                if (lBachTPCClusters < fMinnTPCcls) { AliWarning(Form("Pb / Bach.   track has less than %i TPC clusters ... continue!",fMinnTPCcls)); continue; }
                // - Poor quality related to clusters/findable
                if( lPosTPCFindClusters <= 0 || lNegTPCFindClusters <= 0 || lBachTPCFindClusters ) { AliWarning("Pb / Number of findable cluster <= 0 ... continue!"); continue; }
                if ((lPosTPCClusters/lPosTPCFindClusters)    < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / V0 Pos. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
                if ((lNegTPCClusters/lNegTPCFindClusters)    < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / V0 Neg. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
                if ((lBachTPCClusters/lBachTPCFindClusters)  < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / Bach. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
          }

          //-----------------------------------
          const AliExternalTrackParam *pExtTrack    = pTrackXi->GetInnerParam();
          const AliExternalTrackParam *nExtTrack    = nTrackXi->GetInnerParam();
          const AliExternalTrackParam *bachExtTrack = bachTrackXi->GetInnerParam();
          if (pExtTrack) {
                lInnerWallMomCascDghters[0] = pExtTrack->GetP() * pExtTrack->Charge();
                lTPCSignalCascDghters   [0] = pTrackXi->GetTPCsignal();
          }
          if (nExtTrack) {
                lInnerWallMomCascDghters[1] = nExtTrack->GetP() * nExtTrack->Charge();
                lTPCSignalCascDghters   [1] = nTrackXi->GetTPCsignal();
          }
          if (bachExtTrack) {
                lInnerWallMomCascDghters[2] = bachExtTrack->GetP() * bachExtTrack->Charge();
                lTPCSignalCascDghters   [2] = bachTrackXi->GetTPCsignal();
          }
          etaPos  = pTrackXi->Eta();
          etaNeg  = nTrackXi->Eta();
          etaBach = bachTrackXi->Eta();
          lInvMassLambdaAsCascDghter = xi->GetEffMass(); //This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
          lDcaV0DaughtersXi 	      = xi->GetDcaV0Daughters(); 
          lV0CosineOfPointingAngleXi = xi->GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],
                                                                      lBestPrimaryVtxPos[1],
                                                                      lBestPrimaryVtxPos[2] );
          lDcaV0ToPrimVertexXi = xi->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
          lDcaBachToPrimVertexXi = TMath::Abs( bachTrackXi->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField ) ); //Note: AliExternalTrackParam::GetD returns an algebraic value ...
          xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] ); 
          lV0RadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1] );	
          lDcaPosToPrimVertexXi = TMath::Abs( pTrackXi	->GetD(	lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField ) ); 
          lDcaNegToPrimVertexXi = TMath::Abs( nTrackXi	->GetD(	lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField ) ); 

          //-----------------------------------------
          // - Extra-selection for cascade candidates
          if (fExtraSelections) { //in AliCascadeVertexer
                if (lDcaXiDaughters > 0.3) continue;                                              // in AliCascadeVertexer
                if (lXiCosineOfPointingAngle < 0.999 ) continue;                                  // in AliCascadeVertexer
                if (lDcaV0ToPrimVertexXi < 0.05) continue;                                        // in AliCascadeVertexer
                if (lDcaBachToPrimVertexXi < 0.03) continue;                                      // in AliCascadeVertexer
                if (lDcaV0DaughtersXi > 1.) continue;                                             // in AliV0vertexer
                if ((fCollidingSystem == 0) && (lV0toXiCosineOfPointingAngle < 0.998)) continue; // in AliV0vertexer
                if ((fCollidingSystem == 1) && (lV0CosineOfPointingAngleXi < 0.998)) continue;  // in AliV0vertexer
                if (lDcaPosToPrimVertexXi < 0.1) continue;                                        // in AliV0vertexer
                if (lDcaNegToPrimVertexXi < 0.1) continue;                                        // in AliV0vertexer
                if (lXiRadius < .9) continue;                                                     // in AliCascadeVertexer
                if (lV0RadiusXi < 0.9) continue;                                                  // in AliV0vertexer
          }

          //---------------------------------------------------------------------------------------------------	
          // - Around effective masses. Change mass hypotheses to cover all the possibilities:  Xi-/+, Omega-/+
          if ( bachTrackXi->Charge() < 0 )	{
                //Calculate the effective mass of the Xi- candidate: Xi- hyp. (pdg code 3312)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , 3312); 	
                lInvMassXiMinus = xi->GetEffMassXi();
                //Calculate the effective mass of the Xi- candidate: Omega- hyp. (pdg code 3334)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , 3334); 	
                lInvMassOmegaMinus = xi->GetEffMassXi();
                //Back to "default" hyp. (Xi-)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , 3312); 
          } // end if negative bachelor
          if ( bachTrackXi->Charge() >  0 ) {
                //Calculate the effective mass of the Xi- candidate: Xi+ hyp. (pdg code -3312)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , -3312); 	
                lInvMassXiPlus = xi->GetEffMassXi();
                //Calculate the effective mass of the Xi- candidate: Omega+ hyp. (pdg code -3334)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , -3334); 	
                lInvMassOmegaPlus = xi->GetEffMassXi();
                //Back to "default" hyp. (Xi-)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , -3312); 	
          } // end if positive bachelor

          //--------------------------------
          // - PID on the daughter tracks
          // - Combined PID ->  removed, add when will be used

          // - TPC PID : 3-sigma bands on Bethe-Bloch curve
          //Bachelor
          if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < fTPCPIDsigma) lIsBachelorKaonForTPC = kTRUE;
          if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < fTPCPIDsigma) lIsBachelorPionForTPC = kTRUE;
          //Negative V0 daughter
          if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < fTPCPIDsigma) lIsNegPionForTPC   = kTRUE;
          if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < fTPCPIDsigma) lIsNegProtonForTPC = kTRUE;
          //Positive V0 daughter
          if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < fTPCPIDsigma) lIsPosPionForTPC   = kTRUE;
          if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < fTPCPIDsigma) lIsPosProtonForTPC = kTRUE;
          /*
          const AliExternalTrackParam *pInnerWallTrackXi    = pTrackXi    ->GetInnerParam();
          const AliExternalTrackParam *nInnerWallTrackXi    = nTrackXi    ->GetInnerParam();
          const AliExternalTrackParam *bachInnerWallTrackXi = bachTrackXi ->GetInnerParam();
          if (pInnerWallTrackXi && nInnerWallTrackXi && bachInnerWallTrackXi ) { 
                Double_t pMomInnerWall    = pInnerWallTrackXi   ->GetP();
                Double_t nMomInnerWall    = nInnerWallTrackXi   ->GetP();
                Double_t bachMomInnerWall = bachInnerWallTrackXi->GetP();
                //Bachelor
                if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 3)                              lIsBachelorPionForTPC = kTRUE;
                if (bachMomInnerWall < 0.350  && TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 5) lIsBachelorKaonForTPC = kTRUE;
                if (bachMomInnerWall > 0.350  && TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 3) lIsBachelorKaonForTPC = kTRUE;                
                //Negative V0 daughter
                if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 3  )                           lIsNegPionForTPC   = kTRUE;
                if (nMomInnerWall < 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton ) ) < 5 )   lIsNegProtonForTPC = kTRUE;
                if (nMomInnerWall > 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton ) ) < 3 )   lIsNegProtonForTPC = kTRUE;       
                //Positive V0 daughter
                if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 3 )                            lIsPosPionForTPC   = kTRUE;
                if (pMomInnerWall < 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 5)     lIsPosProtonForTPC = kTRUE;
                if (pMomInnerWall > 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 3)     lIsPosProtonForTPC = kTRUE;
          }*/
		
          //---------------------------------
          // - Extra info for QA (ESD)
          // Miscellaneous pieces of info that may help regarding data quality assessment.
          // Cascade transverse and total momentum
          xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
          lXiTransvMom	= TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY );
          lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY + lXiMomZ*lXiMomZ );
          // V0 total momentum
          xi->GetNPxPyPz(lV0NMomX,lV0NMomY,lV0NMomZ);
          xi->GetPPxPyPz(lV0PMomX,lV0PMomY,lV0PMomZ);
          lV0TotMom = TMath::Sqrt(TMath::Power(lV0NMomX+lV0PMomX,2) + TMath::Power(lV0NMomY+lV0PMomY,2) + TMath::Power(lV0NMomZ+lV0PMomZ,2));
          // Bachelor total momentum
          xi->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
          lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX + lBachMomY*lBachMomY );
          lBachTotMom     = TMath::Sqrt( lBachMomX*lBachMomX + lBachMomY*lBachMomY + lBachMomZ*lBachMomZ );
          lnTrackTransvMom = TMath::Sqrt( lV0NMomX*lV0NMomX   + lV0NMomY*lV0NMomY );
          lpTrackTransvMom = TMath::Sqrt( lV0PMomX*lV0PMomX   + lV0PMomY*lV0PMomY );
          lChargeXi       = xi->Charge();
          lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
          lRapXi    = xi->RapXi();
          lRapOmega = xi->RapOmega();
          lEta      = xi->Eta();
          lTheta    = xi->Theta()*180.0/TMath::Pi();
          lPhi      = xi->Phi()*180.0/TMath::Pi();
          lAlphaXi  = xi->AlphaXi();
          lPtArmXi  = xi->PtArmXi();
	  // Extra-cut = Anti-splitting cut for lambda daughters
          Bool_t kAntiSplittingLambda = kFALSE;	
          if (kAntiSplittingLambda) { // not used
               Double_t lNMomX = 0., lNMomY = 0., lNMomZ = 0.;
               Double_t lPMomX = 0., lPMomY = 0., lPMomZ = 0.;
               xi->GetPPxPyPz(lPMomX, lPMomY, lPMomZ); 
               xi->GetNPxPyPz(lNMomX, lNMomY, lNMomZ); 
               if ( xi->Charge() < 0) {// Xi- or Omega-
                   if (TMath::Abs(lBachTransvMom - TMath::Sqrt( lNMomX*lNMomX + lNMomY*lNMomY )  ) < 0.075) continue;
	       } else {                //Xi+ or Omega+
	           if(TMath::Abs(lBachTransvMom - TMath::Sqrt( lPMomX*lPMomX + lPMomY*lPMomY ) ) < 0.075) continue;
	       }
          }

    } // end of ESD treatment
 
    else if (fAnalysisType == "AOD") {

          // -------------------------------------
          // - Load the cascades from the handler
          const AliAODcascade *xi = lAODevent->GetCascade(iXi);
          if (!xi) continue;
		
          //----------------------------------------------------------------------------        
          // - Assigning the necessary variables for specific AliESDcascade data members		
          lEffMassXi        	   = xi->MassXi(); // default working hypothesis : cascade = Xi- decay
          lDcaXiDaughters	   = xi->DcaXiDaughters();
          lXiCosineOfPointingAngle = xi->CosPointingAngleXi( lBestPrimaryVtxPos[0], 
                                                             lBestPrimaryVtxPos[1], 
		 					     lBestPrimaryVtxPos[2] );
          lPosXi[0] = xi->DecayVertexXiX();
          lPosXi[1] = xi->DecayVertexXiY();
          lPosXi[2] = xi->DecayVertexXiZ();
          lXiRadius = TMath::Sqrt( lPosXi[0]*lPosXi[0] + lPosXi[1]*lPosXi[1] );

          //-------------------------------------------------------------------------------------------------------------------------------
          // - Around the tracks: Bach + V0 (ESD). Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
          AliAODTrack *pTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
          AliAODTrack *nTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
          AliAODTrack *bachTrackXi = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );
          if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
                AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
                continue;
          }
          UInt_t lIdxPosXi  = (UInt_t) TMath::Abs( pTrackXi->GetID() );  
          UInt_t lIdxNegXi  = (UInt_t) TMath::Abs( nTrackXi->GetID() );
          UInt_t lBachIdx   = (UInt_t) TMath::Abs( bachTrackXi->GetID() );
                                                // Care track label can be negative in MC production (linked with the track quality)
                                                // However = normally, not the case for track index ...
          // - Rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
          if (lBachIdx == lIdxNegXi) { AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue; }
          if (lBachIdx == lIdxPosXi) { AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue; }
          // - Get the TPCnumber of cluster for the daughters
          lPosTPCClusters   = pTrackXi->GetTPCClusterInfo(2,1);             // Old for #TPCclusters: ->GetTPCNcls();
          lNegTPCClusters   = nTrackXi->GetTPCClusterInfo(2,1);             // Old for #TPCclusters: ->GetTPCNcls();
          lBachTPCClusters  = bachTrackXi->GetTPCClusterInfo(2,1);          // Old for #TPCclusters: ->GetTPCNcls();
          // - Get the number TPC findable clusters for daughters
          lPosTPCFindClusters   = pTrackXi->GetTPCNclsF();                  // New
          lNegTPCFindClusters   = nTrackXi->GetTPCNclsF();                  // New
          lBachTPCFindClusters  = bachTrackXi->GetTPCNclsF();               // New

          //-------------------------------------
          // - Rejection of a poor quality tracks
          if (fTrackQualityCutTPCrefit) {
                // - Poor quality related to TPCrefit
                if (!(pTrackXi->IsOn(AliAODTrack::kTPCrefit)))    { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
                if (!(nTrackXi->IsOn(AliAODTrack::kTPCrefit)))    { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
                if (!(bachTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
          }
          if (fTrackQualityCutnTPCcls) {
                // - Poor quality related to TPC clusters
                if (lPosTPCClusters  < fMinnTPCcls) { AliWarning(Form("Pb / V0 Pos. track has less than %i TPC clusters ... continue!",fMinnTPCcls)); continue; }
                if (lNegTPCClusters  < fMinnTPCcls) { AliWarning(Form("Pb / V0 Neg. track has less than %i TPC clusters ... continue!",fMinnTPCcls)); continue; }
                if (lBachTPCClusters < fMinnTPCcls) { AliWarning(Form("Pb / Bach. track has less than %i TPC clusters ... continue!",fMinnTPCcls)); continue; }
                // - Poor quality related to clusters/findable
                if( lPosTPCFindClusters <= 0 || lNegTPCFindClusters <= 0 || lBachTPCFindClusters ) { AliWarning("Pb / Number of findable cluster <= 0 ... continue!"); continue; }
                if ((lPosTPCClusters/lPosTPCFindClusters)    < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / V0 Pos. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
                if ((lNegTPCClusters/lNegTPCFindClusters)    < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / V0 Neg. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
                if ((lBachTPCClusters/lBachTPCFindClusters)  < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / Bach. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
          }

          //---------------------------------------
          // - Around the tracks: Bach + V0 (AOD). Necessary variables for AODcascade data members coming from the AODv0 part (inheritance)
          etaPos  = pTrackXi->Eta();
          etaNeg  = nTrackXi->Eta();
          etaBach = bachTrackXi->Eta();
          lChargeXi = xi->ChargeXi();
          if ( lChargeXi < 0) lInvMassLambdaAsCascDghter = xi->MassLambda();
          else 	       lInvMassLambdaAsCascDghter = xi->MassAntiLambda();
          lDcaV0DaughtersXi  	  = xi->DcaV0Daughters(); 
          lDcaV0ToPrimVertexXi   = xi->DcaV0ToPrimVertex();
          lDcaBachToPrimVertexXi = xi->DcaBachToPrimVertex(); 
          lPosV0Xi[0] = xi->DecayVertexV0X();
          lPosV0Xi[1] = xi->DecayVertexV0Y();
          lPosV0Xi[2] = xi->DecayVertexV0Z(); 
          lV0RadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1] );
          lV0CosineOfPointingAngleXi = xi->CosPointingAngle( lBestPrimaryVtxPos ); 
          lDcaPosToPrimVertexXi      = xi->DcaPosToPrimVertex(); 
          lDcaNegToPrimVertexXi      = xi->DcaNegToPrimVertex(); 

          //----------------------------------------------------------------------------------------------------       
          // - Around effective masses. Change mass hypotheses to cover all the possibilities:  Xi-/+, Omega -/+
          if ( lChargeXi < 0 )	lInvMassXiMinus	   = xi->MassXi();
          if ( lChargeXi > 0 )	lInvMassXiPlus 	   = xi->MassXi();
          if ( lChargeXi < 0 )	lInvMassOmegaMinus = xi->MassOmega();
          if ( lChargeXi > 0 )	lInvMassOmegaPlus  = xi->MassOmega();

          //--------------------------------
          // - PID on the daughter tracks
          // - Combined PID ->  removed, add when will be used

          // - TPC PID : 3-sigma bands on Bethe-Bloch curve
          //Bachelor
          if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < fTPCPIDsigma) lIsBachelorKaonForTPC = kTRUE;
          if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < fTPCPIDsigma) lIsBachelorPionForTPC = kTRUE;
          //Negative V0 daughter
          if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < fTPCPIDsigma) lIsNegPionForTPC   = kTRUE;
          if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < fTPCPIDsigma) lIsNegProtonForTPC = kTRUE;
          //Positive V0 daughter
          if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < fTPCPIDsigma) lIsPosPionForTPC   = kTRUE;
          if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < fTPCPIDsigma) lIsPosProtonForTPC = kTRUE;

          //---------------------------------
          // - Extra info for QA (AOD)
          // Miscellaneous pieces of info that may help regarding data quality assessment.
          // Cascade transverse and total momentum	
          lXiMomX = xi->MomXiX();
          lXiMomY = xi->MomXiY();
          lXiMomZ = xi->MomXiZ();
          lXiTransvMom = TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY );
          lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY + lXiMomZ*lXiMomZ );
          Double_t lV0MomX = xi->MomV0X();
          Double_t lV0MomY = xi->MomV0Y();
          Double_t lV0MomZ = xi->MomV0Z();
          lV0TotMom = TMath::Sqrt(TMath::Power(lV0MomX,2)+TMath::Power(lV0MomY,2)+TMath::Power(lV0MomZ,2));
          lBachMomX = xi->MomBachX();
          lBachMomY = xi->MomBachY();
          lBachMomZ = xi->MomBachZ();		
          lBachTransvMom = TMath::Sqrt( lBachMomX*lBachMomX + lBachMomY*lBachMomY );
          lBachTotMom    = TMath::Sqrt( lBachMomX*lBachMomX + lBachMomY*lBachMomY + lBachMomZ*lBachMomZ );
          lV0NMomX = xi->MomNegX();
          lV0NMomY = xi->MomNegY();
          lV0PMomX = xi->MomPosX();
          lV0PMomY = xi->MomPosY();
          lnTrackTransvMom = TMath::Sqrt( lV0NMomX*lV0NMomX   + lV0NMomY*lV0NMomY );
          lpTrackTransvMom = TMath::Sqrt( lV0PMomX*lV0PMomX   + lV0PMomY*lV0PMomY );
          lV0toXiCosineOfPointingAngle = xi->CosPointingAngle( xi->GetDecayVertexXi() );
          lRapXi    = xi->RapXi();
          lRapOmega = xi->RapOmega();
          lEta      = xi->Eta();	              	// Will not work ! need a method Pz(), Py() Px() 
          lTheta    = xi->Theta() *180.0/TMath::Pi();  // in AODcascade.
          lPhi      = xi->Phi()   *180.0/TMath::Pi();  // Here, we will get eta, theta, phi for the V0 ...
          lAlphaXi  = xi->AlphaXi();
          lPtArmXi  = xi->PtArmXi();

    } // end of AOD treatment

    // Cut on pt of the three daughter tracks
    if (lBachTransvMom<fMinPtCutOnDaughterTracks) continue;
    if (lpTrackTransvMom<fMinPtCutOnDaughterTracks) continue;
    if (lnTrackTransvMom<fMinPtCutOnDaughterTracks) continue;
      
    // Cut on pseudorapidity of the three daughter tracks
    if (TMath::Abs(etaBach)>fEtaCutOnDaughterTracks) continue;
    if (TMath::Abs(etaPos)>fEtaCutOnDaughterTracks) continue;
    if (TMath::Abs(etaNeg)>fEtaCutOnDaughterTracks) continue;
      
      
    //----------------------------------
    // Calculate proper lenght for cascade
    Double_t cascadeMass = 0.;
    if ( ( (lChargeXi<0) && lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) ||
         ( (lChargeXi>0) && lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )  ) cascadeMass = 1.321;
    if ( ( (lChargeXi<0) && lIsBachelorKaonForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC ) ||
         ( (lChargeXi>0) && lIsBachelorKaonForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )  ) cascadeMass = 1.672;
    Double_t lctau =  TMath::Sqrt(TMath::Power((lPosXi[0]-lBestPrimaryVtxPos[0]),2)+TMath::Power((lPosXi[1]-lBestPrimaryVtxPos[1]),2)+TMath::Power(( lPosXi[2]-lBestPrimaryVtxPos[2]),2));
    if (lXiTotMom!=0) lctau = lctau*cascadeMass/lXiTotMom;
    else lctau = -1.;

    //-------------------------------------------------
    // Calculate proper lenght for Lambda (reconstructed)
    Float_t lambdaMass = 1.115683; // PDG mass
    Float_t distV0Xi =  TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2)+TMath::Power((lPosV0Xi[2]-lPosXi[2]),2));
    Float_t lctauV0 = -1.;
    if (lV0TotMom!=0) lctauV0 = distV0Xi*lambdaMass/lV0TotMom;
    Float_t distTV0Xi =  TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2));

    //--------------
    /*// For AliEVE      
         if(lChargeXi < 0&& lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) {
             if (lXiTransvMom>2.&&lXiTransvMom<4.&&(lInvMassXiMinus<1.322&&lInvMassXiMinus>1.320)&&(lXiRadius<8.&&lXiRadius>3.)) {
                         // FIXME : Just to know which file is currently open : locate the file containing Xi
                  cout << "Name of the file containing Xi candidate(s) :" 
                       << CurrentFileName() 
                       << " / entry: "     << Entry()
                       << " / in file: "   << lESDevent->GetEventNumberInFile()   // <- Cvetan / From Mihaela: AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetReadEntry();
                       << " AliESDcascade number " << iXi 
                       << " : mass(Xi-) = " << lInvMassXiMinus
                       << " / charge = "   << lChargeXi
                       << " / pt(Casc) = " << lXiTransvMom
                       << " / Decay 2d R(Xi) = " << lXiRadius 
                       << endl;
             }
         }
         if(lChargeXi < 0&& lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) {
             if (lXiTransvMom>2&&lXiTransvMom<4&&(lInvMassOmegaMinus<1.674&&lInvMassOmegaMinus>1.670)&&(lXiRadius<8.&&lXiRadius>3.)) {
                  cout << "Name of the file containing Omega candidate(s) :"
                       << CurrentFileName()
                       << " / entry: "     << Entry()
                       << " / in file: "   << lESDevent->GetEventNumberInFile()   // <- Cvetan / From Mihaela: AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetReadEntry();
                       << " AliESDcascade number " << iXi 
                       << " : mass(Omega-) = " << lInvMassOmegaMinus
                       << " / charge = "   << lChargeXi
                       << " / pt(Casc) = " << lXiTransvMom
                       << " / Decay 2d R(Xi) = " << lXiRadius
                       << endl;

             }
         }*/
          

    // - 
    fHistPosV0TPCClusters->Fill( lPosTPCClusters );
    fHistNegV0TPCClusters->Fill( lNegTPCClusters );
    fHistBachTPCClusters->Fill( lBachTPCClusters );
    f2dHistTPCdEdxOfCascDghters->Fill( lInnerWallMomCascDghters[0] , lTPCSignalCascDghters[0] );
    f2dHistTPCdEdxOfCascDghters->Fill( lInnerWallMomCascDghters[1] , lTPCSignalCascDghters[1] );
    f2dHistTPCdEdxOfCascDghters->Fill( lInnerWallMomCascDghters[2] , lTPCSignalCascDghters[2] );

    //----------------
    //Plot with PID on  
    if ( ( (lChargeXi<0) && lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) ||
         ( (lChargeXi>0) && lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )  ) {
                                       // NOTE : 
                                       // with this condition, it could happen that a cascade candidate satisfies the wrong requirement,
                                       // e.g. one looks at a Xi- candidate for which lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC = kFALSE
                                       // Expectation: it should be excluded, but lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC = kTRUE
                                       // then this bad Xi-candidate will contribute anyway (OR condition).
                                       // Hence: the extra condition on the sign of the Cascade
           //if (TMath::Abs( lInvMassXiMinus-1.3217 ) < 0.010 || TMath::Abs( lInvMassXiPlus-1.3217 ) < 0.010)
           fHistEffMassXi->Fill( lEffMassXi );
           fHistDcaXiDaughters->Fill( lDcaXiDaughters );         	    // Flag CascadeVtxer: Cut Variable e 
           fHistDcaBachToPrimVertex->Fill( lDcaBachToPrimVertexXi   );	    // Flag CascadeVtxer: Cut Variable d
           fHistXiCosineOfPointingAngle->Fill( lXiCosineOfPointingAngle );  // Flag CascadeVtxer: Cut Variable f
           fHistXiRadius->Fill( lXiRadius );	                            // Flag CascadeVtxer: Cut Variable g+h
           fHistMassLambdaAsCascDghter->Fill( lInvMassLambdaAsCascDghter ); // Flag CascadeVtxer: Cut Variable c
           fHistDcaV0DaughtersXi->Fill( lDcaV0DaughtersXi );
           fHistV0CosineOfPointingAngleXi->Fill( lV0CosineOfPointingAngleXi ); 
           fHistV0RadiusXi->Fill( lV0RadiusXi );
           fHistDcaV0ToPrimVertexXi->Fill( lDcaV0ToPrimVertexXi );	    // Flag CascadeVtxer: Cut Variable b
           fHistDcaPosToPrimVertexXi->Fill( lDcaPosToPrimVertexXi );
           fHistDcaNegToPrimVertexXi->Fill( lDcaNegToPrimVertexXi );
           fHistChargeXi->Fill( lChargeXi );
           fHistV0toXiCosineOfPointingAngle->Fill( lV0toXiCosineOfPointingAngle );
           if ( TMath::Abs( lInvMassXiMinus-1.3217 ) < 0.012 || TMath::Abs( lInvMassXiPlus-1.3217 ) < 0.012) { // One InvMass should be different from 0
                fHistXiTransvMom->Fill( lXiTransvMom );
                fHistXiTotMom->Fill( lXiTotMom );
                fHistBachTransvMomXi->Fill( lBachTransvMom );
                fHistBachTotMomXi->Fill( lBachTotMom );
                fHistRapXi->Fill( lRapXi );
                fHistEtaXi->Fill( lEta );
                if (lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) {
                      fHistEtaBachXi->Fill( etaBach );
                      fHistEtaPosXi->Fill( etaPos );
                      fHistEtaNegXi->Fill( etaNeg );
                }
                fHistThetaXi->Fill( lTheta );
                fHistPhiXi->Fill( lPhi );
           }
           if ( TMath::Abs( lInvMassOmegaMinus-1.672 ) < 0.012 || TMath::Abs( lInvMassOmegaPlus-1.672 ) < 0.012 ) { // One InvMass should be different from 0
                fHistRapOmega->Fill( lRapOmega ); 
           }
           f2dHistArmenteros->Fill( lAlphaXi, lPtArmXi );
    } // end with PID ...

    //-----------------------
    // - Invariant mass plots
    //Plots 1D
    if ( lChargeXi < 0 ) {
         fHistMassXiMinus->Fill( lInvMassXiMinus );
         fHistMassOmegaMinus->Fill( lInvMassOmegaMinus );
         f2dHistDcaXiDaughtersvsInvMass->Fill(lDcaXiDaughters,lInvMassXiMinus);
         f2dHistDcaBachToPrimVertexvsInvMass->Fill(lDcaBachToPrimVertexXi,lInvMassXiMinus); 
         f2dHistXiCosineOfPointingAnglevsInvMass->Fill(lXiCosineOfPointingAngle,lInvMassXiMinus);
         f2dHistMassLambdaAsCascDghtervsInvMass->Fill(lInvMassLambdaAsCascDghter,lInvMassXiMinus);
         f2dHistDcaV0DaughtersXivsInvMass->Fill(lDcaV0DaughtersXi,lInvMassXiMinus);
         f2dHistDcaV0ToPrimVertexXivsInvMass->Fill(lDcaV0ToPrimVertexXi,lInvMassXiMinus);
    }
    if ( lChargeXi > 0 ) {
      fHistMassXiPlus->Fill( lInvMassXiPlus );
      fHistMassOmegaPlus->Fill( lInvMassOmegaPlus );
    }
    //Plots 2D, 3D
    if ( lChargeXi < 0 ) {
      f2dHistEffMassLambdaVsEffMassXiMinus->Fill( lInvMassLambdaAsCascDghter, lInvMassXiMinus ); 
      f2dHistEffMassXiVsEffMassOmegaMinus ->Fill( lInvMassXiMinus, lInvMassOmegaMinus );
      f2dHistXiRadiusVsEffMassXiMinus     ->Fill( lXiRadius, lInvMassXiMinus );
      f2dHistXiRadiusVsEffMassOmegaMinus  ->Fill( lXiRadius, lInvMassOmegaMinus );
    } else {
      f2dHistEffMassLambdaVsEffMassXiPlus ->Fill( lInvMassLambdaAsCascDghter, lInvMassXiPlus );
      f2dHistEffMassXiVsEffMassOmegaPlus  ->Fill( lInvMassXiPlus, lInvMassOmegaPlus );
      f2dHistXiRadiusVsEffMassXiPlus      ->Fill( lXiRadius, lInvMassXiPlus);
      f2dHistXiRadiusVsEffMassOmegaPlus   ->Fill( lXiRadius, lInvMassOmegaPlus );
    }

    //---------------------------------------------	
    // - Filling the AliCFContainers related to PID
    Double_t lContainerPIDVars[3] = {0.0};
    // Xi Minus		
    if ( lChargeXi < 0 ) {
          lContainerPIDVars[0] = lXiTransvMom;
          lContainerPIDVars[1] = lInvMassXiMinus;
          lContainerPIDVars[2] = lRapXi;
          //No PID
          fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 0); // No PID
	  //TPC PID
          if ( lIsBachelorPionForTPC )                                           fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track		
          if ( lIsBachelorPionForTPC && lIsPosProtonForTPC )                     fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
          if ( lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	  //Combined PID
          if ( lIsBachelorPion )                                      fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor	
          if ( lIsBachelorPion && lIsPosInXiProton )                  fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
          if (lIsBachelorPion && lIsPosInXiProton && lIsNegInXiPion ) fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
    lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; 	
    // Xi Plus		
    if ( lChargeXi > 0 ) {
          lContainerPIDVars[0] = lXiTransvMom;
          lContainerPIDVars[1] = lInvMassXiPlus;
          lContainerPIDVars[2] = lRapXi;
	  //No PID
          fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 0); // No PID
          //TPC PID
          if ( lIsBachelorPionForTPC )                                           fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
          if ( lIsBachelorPionForTPC && lIsNegProtonForTPC )                     fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
          if ( lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC ) fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	  //Combined PID
          if ( lIsBachelorPion )                                      fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
          if ( lIsBachelorPion && lIsNegInXiProton )                  fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
          if (lIsBachelorPion && lIsNegInXiProton && lIsPosInXiPion ) fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
    lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; 
    // Omega Minus		
    if ( lChargeXi < 0 ) {
          lContainerPIDVars[0] = lXiTransvMom;
          lContainerPIDVars[1] = lInvMassOmegaMinus;
          lContainerPIDVars[2] = lRapOmega;
	  //No PID
          fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 0); // No PID
       	  //TPC PID
          if ( lIsBachelorKaonForTPC )                                           fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
          if ( lIsBachelorKaonForTPC && lIsPosProtonForTPC )                     fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
          if ( lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	  //Combined PID
          if ( lIsBachelorKaon )                                            fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
          if ( lIsBachelorKaon && lIsPosInOmegaProton )     	            fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
          if (lIsBachelorKaon && lIsPosInOmegaProton && lIsNegInOmegaPion ) fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
    lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.;
    // Omega Plus		
    if ( lChargeXi > 0 ) {
      lContainerPIDVars[0] = lXiTransvMom;
      lContainerPIDVars[1] = lInvMassOmegaPlus;
      lContainerPIDVars[2] = lRapOmega;
	// No PID
      fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 0); // No PID
	// TPC PID
      if ( lIsBachelorKaonForTPC  )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
      if( lIsBachelorKaonForTPC && 
	  lIsNegProtonForTPC     )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
      if ( lIsBachelorKaonForTPC && 
	   lIsNegProtonForTPC    && 
	   lIsPosPionForTPC       )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	// Combined PID
      if ( lIsBachelorKaon        )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
      if ( lIsBachelorKaon       && 
           lIsNegInOmegaProton    )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
      if (lIsBachelorKaon     && 
	  lIsNegInOmegaProton && 
	  lIsPosInOmegaPion    )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
        	
    //--------------------------------------------------------------------
    // Filling the AliCFContainer (optimisation of topological selections)
    Double_t lContainerCutVars[20] = {0.0};
                        
    lContainerCutVars[0]  = lDcaXiDaughters;
    lContainerCutVars[1]  = lDcaBachToPrimVertexXi;
    lContainerCutVars[2]  = lXiCosineOfPointingAngle;
    lContainerCutVars[3]  = lXiRadius;
    lContainerCutVars[4]  = lInvMassLambdaAsCascDghter;
    lContainerCutVars[5]  = lDcaV0DaughtersXi;
    if      (fCollidingSystem == 0) lContainerCutVars[6]  = lV0toXiCosineOfPointingAngle;
    else if (fCollidingSystem == 1) lContainerCutVars[6]  = lV0CosineOfPointingAngleXi;
    lContainerCutVars[7]  = lV0RadiusXi;
    lContainerCutVars[8]  = lDcaV0ToPrimVertexXi;	
    lContainerCutVars[9]  = lDcaPosToPrimVertexXi;
    lContainerCutVars[10] = lDcaNegToPrimVertexXi;
    lContainerCutVars[13] = lXiTransvMom;
    lContainerCutVars[16] = lctau;
    lContainerCutVars[17] = lctauV0;
    lContainerCutVars[18] = distTV0Xi;
    if ( lChargeXi < 0 ) {
         lContainerCutVars[11] = lInvMassXiMinus;
         lContainerCutVars[12] = lInvMassOmegaMinus;
         lContainerCutVars[14] = lRapXi;
         lContainerCutVars[15] = -1.;
         lContainerCutVars[19] = TMath::Abs(lInnerWallMomCascDghters[0]);
         if (lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) fCFContCascadeCuts->Fill(lContainerCutVars,0); // for Xi-
         lContainerCutVars[11] = lInvMassXiMinus;
         lContainerCutVars[12] = lInvMassOmegaMinus;
         lContainerCutVars[14] = -1.;
         lContainerCutVars[15] = lRapOmega;
         lContainerCutVars[19] = TMath::Abs(lInnerWallMomCascDghters[0]); 
         if (lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) fCFContCascadeCuts->Fill(lContainerCutVars,2); // for Omega-
    } else {
         lContainerCutVars[11] = lInvMassXiPlus;
         lContainerCutVars[12] = lInvMassOmegaPlus; 
         lContainerCutVars[14] = lRapXi;
         lContainerCutVars[15] = -1.; 
         lContainerCutVars[19] = TMath::Abs(lInnerWallMomCascDghters[1]);
         if (lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC) fCFContCascadeCuts->Fill(lContainerCutVars,1); // for Xi+
         lContainerCutVars[11] = lInvMassXiPlus;
         lContainerCutVars[12] = lInvMassOmegaPlus;
         lContainerCutVars[14] = -1.;
         lContainerCutVars[15] = lRapOmega;
         lContainerCutVars[19] = TMath::Abs(lInnerWallMomCascDghters[1]);
         if (lIsBachelorKaonForTPC && lIsNegProtonForTPC && lIsPosPionForTPC) fCFContCascadeCuts->Fill(lContainerCutVars,3); // for Omega+ 
    }                 
  } //end of the Cascade loop (ESD or AOD)
    
  // Post output data.
 PostData(1, fListHistCascade);
 PostData(2, fCFContCascadePIDXiMinus);
 PostData(3, fCFContCascadePIDXiPlus);
 PostData(4, fCFContCascadePIDOmegaMinus);
 PostData(5, fCFContCascadePIDOmegaPlus);
 PostData(6, fCFContCascadeCuts);
}

//________________________________________________________________________
Int_t AliAnalysisTaskCheckCascadepp::DoESDTrackWithTPCrefitMultiplicity(const AliESDEvent *lESDevent) {
    // Checking the number of tracks with TPCrefit for each event
    // Needed for a rough assessment of the event multiplicity
        
    Int_t nTrackWithTPCrefitMultiplicity = 0;
    for (Int_t iTrackIdx = 0; iTrackIdx < (InputEvent())->GetNumberOfTracks(); iTrackIdx++) {
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

    } // end loop over all event tracks
    return  nTrackWithTPCrefitMultiplicity;
}


//________________________________________________________________________
void AliAnalysisTaskCheckCascadepp::Terminate(Option_t *) 
{
  
  
}
