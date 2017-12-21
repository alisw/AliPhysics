/***************************************************************          *
 *  Authors : Antonin Maire, Boris Hippolyte 
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
//	      AliAnalysisTaskCheckPerformanceCascadePbPb class
//            This task is for a performance study of cascade identification in PbPb.
//            It works with MC info and ESD/AOD.
//            Origin   : AliAnalysisTaskCheckPerformanceCascade class by A. Maire Nov2010, antonin.maire@ires.in2p3.fr
//            Modified for PbPb analysis: M. Nicassio Feb2011, maria.nicassio@ba.infn.it:
//                        - physics selection moved to the run.C macro
//                        - added centrality selection and possibility to select events in nTracks ranges 
//                        - added new histograms 
//                        - modified binning of some histograms and containers 
//                        - flag to enable CF container usage 
//                        - check in the destructor for CAF usage
//                        - flag for acceptance cut in the MC part
//                        - in the MC particle selection IsPhysicalPrimary added and number of particles taken as appropriate for HIJING 
//                          (however for cascades one gets the same if runs on Nprimaries in the stack and does not check IsPhysicalPrimary)
//                        - automatic settings for PID 
//                        - selection of injected cascades and HIJING cascades (kind of "bug" in method IsFromBGEvent())
//                        - added proper time histograms for cascades and lambdas
//                        - cos of PA V0 wrt Xi vertex and not primary vertex  
//                        - distance xi-V0 added in the container
//                        - AOD analysis developed (January 2012)
//                   
//            Adapted to pp 2.76 analysis: D. Colella, domenico.colella@ba.infn.it (Nov. 2012):
//                        - added new and removed other histograms 
//                        - Physics selection moved here (mainly for normalization in the efficiency calcuation)
//                        - Centrality selection deleted
//                        - 3DHisto denominator moved before any event selection for Normalization
//                        - injected and natural part of MC selection removed
//
//            Adapted to pPb 5.02 TeV analysis: D. Colella, domenico.colella@cern.ch (Sep. 2014)
//                        - Added the parameter fCollidingSystem, to distingish between pp and pPb procedures
//
//            Adapted to pp 13 TeV analysis: D. Colella, domenico.colella@cern.ch (Aug. 2015)
//                        - Clarify the usage of SDD selection after the introduction of fCollidingSystem
//
//            Improvement for pp 13 TeV analysis: D. Colella, domenico.colella@cern.ch
//               Feb 2016
//               - Added the possibility to select events in the class kINT7 for the pp@13TeV analysis through "fkSwitchINT7"
//               - Revision of the event selections for the pp@13TeV analysis
//               March 2016
//               - Checks on the normalization and the efficiency denominator added
//               May 2016
//               - move the "fCollidingSystem" from TString to Int_t
//               - add fSuffix for the LEGO train usage
//               - More checks on the event selections
//
//-----------------------------------------------------------------


#include <Riostream.h>

#include "TList.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TMath.h"
#include "AliLog.h"
#include "AliHeader.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliCFContainer.h"
#include "AliESDVZERO.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDEvent.h"
#include "AliESDcascade.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h" 
#include "AliAnalysisTaskCheckPerformanceCascadepp.h"
#include "AliAnalysisUtils.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskCheckPerformanceCascadepp)



//________________________________________________________________________________________
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// - DEFAULT CONSTRUCTOR
AliAnalysisTaskCheckPerformanceCascadepp::AliAnalysisTaskCheckPerformanceCascadepp() 
: AliAnalysisTaskSE(), // <- take care to AliAnalysisTask( empty )
    fAnalysisType                   ("ESD"),
    fCollidingSystem                (0),
    fkTriggerClass                  (AliVEvent::kINT7),
    fApplyEvSelSDDstatus            (kTRUE),
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
    fkExtraSelections               (0),
    fVtxRangeMax                    (0),
    fVtxRangeMin                    (0),
    fApplyAccCut                    (0),
    fMinPtCutOnDaughterTracks       (0),
    fEtaCutOnDaughterTracks         (0),
    fSPDPileUpminContributors       (3),
    fTPCPIDsigma                    (4),
    fSuffix                         (""),

    // - Plots initialisation
    // - List
    fListHistCascade(0),
    // - General Plots
      // -- Cascade multiplicity plots
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
      // -- Cascade multiplicity per event plots
      fHistnXiPlusPerEvTot(0),                  // After any event selections, in all the eta and pt range
      fHistnXiMinusPerEvTot(0),                 // After any event selections, in all the eta and pt range
      fHistnOmegaPlusPerEvTot(0),               // After any event selections, in all the eta and pt range
      fHistnOmegaMinusPerEvTot(0),              // After any event selections, in all the eta and pt range
      fHistnXiPlusPerEv(0),                     // After any event selections, if (fApplyAccCut) in the detector acceptance and over a pt minimum
      fHistnXiMinusPerEv(0),                    // After any event selections, if (fApplyAccCut) in the detector acceptance and over a pt minimum
      fHistnOmegaPlusPerEv(0),                  // After any event selections, if (fApplyAccCut) in the detector acceptance and over a pt minimum
      fHistnOmegaMinusPerEv(0),                 // After any event selections, if (fApplyAccCut) in the detector acceptance and over a pt minimum
      fHistnAssoXiMinus(0),                     // For the Reconstructed-Associated cascades 
      fHistnAssoXiPlus(0),                      // For the Reconstructed-Associated cascades 
      fHistnAssoOmegaMinus(0),                  // For the Reconstructed-Associated cascades 
      fHistnAssoOmegaPlus(0),                   // For the Reconstructed-Associated cascades 
      // -- Vertex position plots (BestVertex)
      fHistPVx(0),                              // After any selections but before |Z| < 10 cm
      fHistPVy(0),                              // After any selections but before |Z| < 10 cm
      fHistPVz(0),                              // After any selections but before |Z| < 10 cm
      fHistPVxAnalysis(0),                      // After any event selections
      fHistPVyAnalysis(0),                      // After any event selections
      fHistPVzAnalysis(0),                      // After any event selections
      // TPC cluster distributions for daughters
      fHistPosV0TPCClusters(0),
      fHistNegV0TPCClusters(0),
      fHistBachTPCClusters(0),
      // -- Plots needed for efficiency denominator calculation
        // - Step A) filled before all the selections 
      f3dHistGenPtVsGenYvsNtracksXiMinus_A(0),     
      f3dHistGenPtVsGenYvsNtracksXiPlus_A(0),      
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_A(0),  
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_A(0),  
        // - Step B) filled after the pre-trigger selection (DAQ incomplete, SPD background, Pile-up) --> Needed for the efficiency calculation in method 1
      f3dHistGenPtVsGenYvsNtracksXiMinus_B(0),     
      f3dHistGenPtVsGenctauvsYXiMinus_B(0),       
      f3dHistGenPtVsGenYvsNtracksXiPlus_B(0),     
      f3dHistGenPtVsGenctauvsYXiPlus_B(0),        
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_B(0), 
      f3dHistGenPtVsGenctauvsYOmegaMinus_B(0),    
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_B(0),  
      f3dHistGenPtVsGenctauvsYOmegaPlus_B(0),
        // - Step C) filled after the trigger selection (Physics selection)
      f3dHistGenPtVsGenYvsNtracksXiMinus_C(0),    
      f3dHistGenPtVsGenYvsNtracksXiPlus_C(0),     
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_C(0),                           
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_C(0),       
        // - Step D) filled after the revertexing of the V0 and cascades
      f3dHistGenPtVsGenYvsNtracksXiMinus_D(0),    
      f3dHistGenPtVsGenYvsNtracksXiPlus_D(0),     
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_D(0),  
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_D(0),   
        // - Step E) filled after the request of the presence of the vertex SPD and not only the one from the tracks
      f3dHistGenPtVsGenYvsNtracksXiMinus_E(0),     
      f3dHistGenPtVsGenYvsNtracksXiPlus_E(0),     
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_E(0),  
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_E(0), 
        // - Step F) filled after the request on the vertex resolution and dispersion
      f3dHistGenPtVsGenYvsNtracksXiMinus_F(0),                                 
      f3dHistGenPtVsGenYvsNtracksXiPlus_F(0),     
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_F(0),                              
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_F(0),   
        // - Step G) filled after the request on the vertices proximity
      f3dHistGenPtVsGenYvsNtracksXiMinus_G(0),                                 
      f3dHistGenPtVsGenYvsNtracksXiPlus_G(0),     
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_G(0),                              
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_G(0),   
        // - Step H) filled after the request on |Zpv| < 10 cm that means after all the event selections
      f3dHistGenPtVsGenYvsNtracksXiMinus_H(0),                                                                                                
      f3dHistGenPtVsGenctauvsYXiMinus_H(0),       
      f3dHistGenPtVsGenYvsNtracksXiPlus_H(0),     
      f3dHistGenPtVsGenctauvsYXiPlus_H(0),        
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_H(0), 
      f3dHistGenPtVsGenctauvsYOmegaMinus_H(0),    
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_H(0),  
      f3dHistGenPtVsGenctauvsYOmegaPlus_H(0),  
      // - Plots for generated cascades (after all the event selections)
      // -- xi minus
      fHistEtaGenCascXiMinus(0),                // In all the eta and pt range (as they are generated)
      fHistThetaGenCascXiMinus(0),              // In all the eta and pt range (as they are generated)
      f2dHistGenPtVsGenYFdblXiMinus(0),         // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaLambdaXiMinus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBachXiMinus(0),                 // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaMesDghterXiMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBarDghterXiMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaLambdaXiMinus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBachXiMinus(0),                 // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaMesDghterXiMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBarDghterXiMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBachXiMinus(0),                    // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtMesDghterXiMinus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBarDghterXiMinus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      // -- xi plus
      fHistEtaGenCascXiPlus(0),                 // In all the eta and pt range (as they are generated)
      fHistThetaGenCascXiPlus(0),               // In all the eta and pt range (as they are generated)
      f2dHistGenPtVsGenYFdblXiPlus(0),          // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaLambdaXiPlus(0),                // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBachXiPlus(0),                  // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaMesDghterXiPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBarDghterXiPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaLambdaXiPlus(0),                // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBachXiPlus(0),                  // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaMesDghterXiPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBarDghterXiPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBachXiPlus(0),                     // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtMesDghterXiPlus(0),                // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBarDghterXiPlus(0),                // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      // -- omega minus
      fHistEtaGenCascOmegaMinus(0),             // In all the eta and pt range (as they are generated)
      fHistThetaGenCascOmegaMinus(0),           // In all the eta and pt range (as they are generated)
      f2dHistGenPtVsGenYFdblOmegaMinus(0),      // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaLambdaOmegaMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBachOmegaMinus(0),              // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaMesDghterOmegaMinus(0),         // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBarDghterOmegaMinus(0),         // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaLambdaOmegaMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBachOmegaMinus(0),              // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaMesDghterOmegaMinus(0),         // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBarDghterOmegaMinus(0),         // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBachOmegaMinus(0),                 // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtMesDghterOmegaMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBarDghterOmegaMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      // -- omega plus
      fHistEtaGenCascOmegaPlus(0),              // In all the eta and pt range (as they are generated)
      fHistThetaGenCascOmegaPlus(0),            // In all the eta and pt range (as they are generated)
      f2dHistGenPtVsGenYFdblOmegaPlus(0),       // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaLambdaOmegaPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBachOmegaPlus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaMesDghterOmegaPlus(0),          // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBarDghterOmegaPlus(0),          // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaLambdaOmegaPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBachOmegaPlus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaMesDghterOmegaPlus(0),          // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBarDghterOmegaPlus(0),          // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBachOmegaPlus(0),                  // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtMesDghterOmegaPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBarDghterOmegaPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
    // - Plots for associated to MC cascade (after all the event selections)
      // -- PID Probability versus MC Pt(bachelor track)
      f2dHistPIDprobaKaonVsMCPtBach(0), 
      f2dHistPIDprobaPionVsMCPtBach(0),
      // -- Invariant mass distribution for the cascade candidates associated with MC
      fHistAsMCMassXiMinus(0),		
      fHistAsMCMassXiPlus(0),		
      fHistAsMCMassOmegaMinus(0),
      fHistAsMCMassOmegaPlus(0),
      // -- Generated Pt Vs generated y, for the cascade candidates associated with MC + Info Comb. PID
      f2dHistAsMCandCombPIDGenPtVsGenYXiMinus(0),
      f2dHistAsMCandCombPIDGenPtVsGenYXiPlus(0),
      f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus(0),
      f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus(0),
      // -- Generated Pt Vs generated y, for the cascade candidates associated with MC
      f2dHistAsMCGenPtVsGenYXiMinus(0),
      f2dHistAsMCGenPtVsGenYXiPlus(0),
      f2dHistAsMCGenPtVsGenYOmegaMinus(0),
      f2dHistAsMCGenPtVsGenYOmegaPlus(0),
      // -- Generated Eta of the cascade candidates associated with MC
      fHistAsMCGenEtaXiMinus(0),
      fHistAsMCGenEtaXiPlus(0),
      fHistAsMCGenEtaOmegaMinus(0),
      fHistAsMCGenEtaOmegaPlus(0),
      // -- Resolution in Pt as function of generated Pt
      f2dHistAsMCResPtXiMinus(0),		
      f2dHistAsMCResPtXiPlus(0),		
      f2dHistAsMCResPtOmegaMinus(0),
      f2dHistAsMCResPtOmegaPlus(0),	
      // -- Resolution in R(2D) as function of generated R
      f2dHistAsMCResRXiMinus(0),		
      f2dHistAsMCResRXiPlus(0),		
      f2dHistAsMCResROmegaMinus(0),
      f2dHistAsMCResROmegaPlus(0),
      // -- Resolution in phi as function of generated Pt
      f2dHistAsMCResPhiXiMinus(0),
      f2dHistAsMCResPhiXiPlus(0),
      f2dHistAsMCResPhiOmegaMinus(0),
      f2dHistAsMCResPhiOmegaPlus(0),
      // -- Correlation between proton (antiproton) daughter MC pt and Xi/Omega MC pt (to apply Geant/Fluka correction)
      f2dHistAsMCptProtonMCptXiMinus(0),
      f2dHistAsMCptAntiprotonMCptXiPlus(0),
      f2dHistAsMCptProtonMCptOmegaMinus(0),
      f2dHistAsMCptAntiprotonMCptOmegaPlus(0),
      // -- Containers                       
      fCFContCascadePIDAsXiMinus(0),
      fCFContCascadePIDAsXiPlus(0),
      fCFContCascadePIDAsOmegaMinus(0),
      fCFContCascadePIDAsOmegaPlus(0),
      fCFContAsCascadeCuts(0)

    //____Dummy costructor____
    {
        for(Int_t iV0selIdx   = 0; iV0selIdx   < 7; iV0selIdx++   ) { fV0Sels          [iV0selIdx   ] = -1.; }
        for(Int_t iCascSelIdx = 0; iCascSelIdx < 8; iCascSelIdx++ ) { fCascSels        [iCascSelIdx ] = -1.; }
    }
     
        
     
//____________________________________________________________________________________________
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// - NON DEFAULT CONSTRUCTOR
AliAnalysisTaskCheckPerformanceCascadepp::AliAnalysisTaskCheckPerformanceCascadepp(const char *name) 
  : AliAnalysisTaskSE(name),
    fAnalysisType                   ("ESD"),
    fCollidingSystem                (0),
    fkTriggerClass                  (AliVEvent::kINT7),
    fApplyEvSelSDDstatus            (kTRUE),
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
    fkExtraSelections               (0),
    fVtxRangeMax                    (0),
    fVtxRangeMin                    (0),
    fApplyAccCut                    (0),
    fMinPtCutOnDaughterTracks       (0),
    fEtaCutOnDaughterTracks         (0),
    fSPDPileUpminContributors       (3),
    fTPCPIDsigma                    (4),
    fSuffix                         (""),

    // - Plots initialisation
    fListHistCascade(0),
    // - General Plots
      // -- Cascade multiplicity plots
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
      // -- Cascade multiplicity per event plots
      fHistnXiPlusPerEvTot(0),                  // After any event selections, in all the eta and pt range
      fHistnXiMinusPerEvTot(0),                 // After any event selections, in all the eta and pt range
      fHistnOmegaPlusPerEvTot(0),               // After any event selections, in all the eta and pt range
      fHistnOmegaMinusPerEvTot(0),              // After any event selections, in all the eta and pt range
      fHistnXiPlusPerEv(0),                     // After any event selections, if (fApplyAccCut) in the detector acceptance and over a pt minimum
      fHistnXiMinusPerEv(0),                    // After any event selections, if (fApplyAccCut) in the detector acceptance and over a pt minimum
      fHistnOmegaPlusPerEv(0),                  // After any event selections, if (fApplyAccCut) in the detector acceptance and over a pt minimum
      fHistnOmegaMinusPerEv(0),                 // After any event selections, if (fApplyAccCut) in the detector acceptance and over a pt minimum
      fHistnAssoXiMinus(0),                     // For the Reconstructed-Associated cascades 
      fHistnAssoXiPlus(0),                      // For the Reconstructed-Associated cascades 
      fHistnAssoOmegaMinus(0),                  // For the Reconstructed-Associated cascades 
      fHistnAssoOmegaPlus(0),                   // For the Reconstructed-Associated cascades 
      // -- Vertex position plots (BestVertex)
      fHistPVx(0),                              // After any selections but before |Z| < 10 cm
      fHistPVy(0),                              // After any selections but before |Z| < 10 cm
      fHistPVz(0),                              // After any selections but before |Z| < 10 cm
      fHistPVxAnalysis(0),                      // After any event selections
      fHistPVyAnalysis(0),                      // After any event selections
      fHistPVzAnalysis(0),                      // After any event selections
      // TPC cluster distributions for daughters
      fHistPosV0TPCClusters(0),
      fHistNegV0TPCClusters(0),
      fHistBachTPCClusters(0),
      // -- Plots needed for efficiency denominator calculation
        // - Step A) filled before all the selections 
      f3dHistGenPtVsGenYvsNtracksXiMinus_A(0),      
      f3dHistGenPtVsGenYvsNtracksXiPlus_A(0),                                  
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_A(0),                              
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_A(0),                               
        // - Step B) filled after the pre-trigger selection (DAQ incomplete, SPD background, Pile-up) --> Needed for the efficiency calculation in method 1
      f3dHistGenPtVsGenYvsNtracksXiMinus_B(0),                                                                                                
      f3dHistGenPtVsGenctauvsYXiMinus_B(0),       
      f3dHistGenPtVsGenYvsNtracksXiPlus_B(0),     
      f3dHistGenPtVsGenctauvsYXiPlus_B(0),        
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_B(0), 
      f3dHistGenPtVsGenctauvsYOmegaMinus_B(0),    
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_B(0),  
      f3dHistGenPtVsGenctauvsYOmegaPlus_B(0),     
        // - Step C) filled after the trigger selection (Physics selection)
      f3dHistGenPtVsGenYvsNtracksXiMinus_C(0),
      f3dHistGenPtVsGenYvsNtracksXiPlus_C(0),
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_C(0),
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_C(0),  
        // - Step D) filled after the revertexing of the V0 and cascades
      f3dHistGenPtVsGenYvsNtracksXiMinus_D(0),    
      f3dHistGenPtVsGenYvsNtracksXiPlus_D(0),     
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_D(0),                           
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_D(0),                            
        // - Step E) filled after the request of the presence of the vertex SPD and not only the one from the tracks
      f3dHistGenPtVsGenYvsNtracksXiMinus_E(0),                                 
      f3dHistGenPtVsGenYvsNtracksXiPlus_E(0),     
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_E(0),                              
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_E(0),    
        // - Step F) filled after the request on the vertex resolution and dispersion
      f3dHistGenPtVsGenYvsNtracksXiMinus_F(0),
      f3dHistGenPtVsGenYvsNtracksXiPlus_F(0),
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_F(0),
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_F(0),   
        // - Step G) filled after the request on the vertices proximity
      f3dHistGenPtVsGenYvsNtracksXiMinus_G(0),
      f3dHistGenPtVsGenYvsNtracksXiPlus_G(0),
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_G(0),
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_G(0),   
        // - Step H) filled after the request on |Zpv| < 10 cm that means after all the event selections
      f3dHistGenPtVsGenYvsNtracksXiMinus_H(0),
      f3dHistGenPtVsGenctauvsYXiMinus_H(0),
      f3dHistGenPtVsGenYvsNtracksXiPlus_H(0),
      f3dHistGenPtVsGenctauvsYXiPlus_H(0),
      f3dHistGenPtVsGenYvsNtracksOmegaMinus_H(0),
      f3dHistGenPtVsGenctauvsYOmegaMinus_H(0),
      f3dHistGenPtVsGenYvsNtracksOmegaPlus_H(0),
      f3dHistGenPtVsGenctauvsYOmegaPlus_H(0),
      // - Plots for generated cascades (after all the event selections)
      // -- xi minus
      fHistEtaGenCascXiMinus(0),                // In all the eta and pt range (as they are generated)
      fHistThetaGenCascXiMinus(0),              // In all the eta and pt range (as they are generated)
      f2dHistGenPtVsGenYFdblXiMinus(0),         // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaLambdaXiMinus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBachXiMinus(0),                 // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaMesDghterXiMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBarDghterXiMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaLambdaXiMinus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBachXiMinus(0),                 // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaMesDghterXiMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBarDghterXiMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBachXiMinus(0),                    // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtMesDghterXiMinus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBarDghterXiMinus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      // -- xi plus
      fHistEtaGenCascXiPlus(0),                 // In all the eta and pt range (as they are generated)
      fHistThetaGenCascXiPlus(0),               // In all the eta and pt range (as they are generated)
      f2dHistGenPtVsGenYFdblXiPlus(0),          // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaLambdaXiPlus(0),                // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBachXiPlus(0),                  // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaMesDghterXiPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBarDghterXiPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaLambdaXiPlus(0),                // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBachXiPlus(0),                  // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaMesDghterXiPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBarDghterXiPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBachXiPlus(0),                     // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtMesDghterXiPlus(0),                // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBarDghterXiPlus(0),                // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      // -- omega minus
      fHistEtaGenCascOmegaMinus(0),             // In all the eta and pt range (as they are generated)
      fHistThetaGenCascOmegaMinus(0),           // In all the eta and pt range (as they are generated)
      f2dHistGenPtVsGenYFdblOmegaMinus(0),      // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaLambdaOmegaMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBachOmegaMinus(0),              // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaMesDghterOmegaMinus(0),         // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBarDghterOmegaMinus(0),         // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaLambdaOmegaMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBachOmegaMinus(0),              // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaMesDghterOmegaMinus(0),         // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBarDghterOmegaMinus(0),         // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBachOmegaMinus(0),                 // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtMesDghterOmegaMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBarDghterOmegaMinus(0),            // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      // -- omega plus
      fHistEtaGenCascOmegaPlus(0),              // In all the eta and pt range (as they are generated)
      fHistThetaGenCascOmegaPlus(0),            // In all the eta and pt range (as they are generated)
      f2dHistGenPtVsGenYFdblOmegaPlus(0),       // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaLambdaOmegaPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBachOmegaPlus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaMesDghterOmegaPlus(0),          // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistThetaBarDghterOmegaPlus(0),          // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaLambdaOmegaPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBachOmegaPlus(0),               // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaMesDghterOmegaPlus(0),          // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistEtaBarDghterOmegaPlus(0),          // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBachOmegaPlus(0),                  // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtMesDghterOmegaPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
      fHistPtBarDghterOmegaPlus(0),             // if (fApplyAccCut) In the detector acceptance and over a pt minimum (Findable particle)
    // - Plots for associated to MC cascade (after all the event selections)
      // -- PID Probability versus MC Pt(bachelor track)
      f2dHistPIDprobaKaonVsMCPtBach(0),
      f2dHistPIDprobaPionVsMCPtBach(0),
      // -- Invariant mass distribution for the cascade candidates associated with MC
      fHistAsMCMassXiMinus(0),
      fHistAsMCMassXiPlus(0),
      fHistAsMCMassOmegaMinus(0),
      fHistAsMCMassOmegaPlus(0),
      // -- Generated Pt Vs generated y, for the cascade candidates associated with MC + Info Comb. PID
      f2dHistAsMCandCombPIDGenPtVsGenYXiMinus(0),
      f2dHistAsMCandCombPIDGenPtVsGenYXiPlus(0),
      f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus(0),
      f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus(0),
      // -- Generated Pt Vs generated y, for the cascade candidates associated with MC
      f2dHistAsMCGenPtVsGenYXiMinus(0),
      f2dHistAsMCGenPtVsGenYXiPlus(0),
      f2dHistAsMCGenPtVsGenYOmegaMinus(0),
      f2dHistAsMCGenPtVsGenYOmegaPlus(0),
      // -- Generated Eta of the cascade candidates associated with MC
      fHistAsMCGenEtaXiMinus(0),
      fHistAsMCGenEtaXiPlus(0),
      fHistAsMCGenEtaOmegaMinus(0),
      fHistAsMCGenEtaOmegaPlus(0),
      // -- Resolution in Pt as function of generated Pt
      f2dHistAsMCResPtXiMinus(0),
      f2dHistAsMCResPtXiPlus(0),
      f2dHistAsMCResPtOmegaMinus(0),
      f2dHistAsMCResPtOmegaPlus(0),
      // -- Resolution in R(2D) as function of generated R
      f2dHistAsMCResRXiMinus(0),
      f2dHistAsMCResRXiPlus(0),
      f2dHistAsMCResROmegaMinus(0),
      f2dHistAsMCResROmegaPlus(0),
      // -- Resolution in phi as function of generated Pt
      f2dHistAsMCResPhiXiMinus(0),
      f2dHistAsMCResPhiXiPlus(0),
      f2dHistAsMCResPhiOmegaMinus(0),
      f2dHistAsMCResPhiOmegaPlus(0),
      // -- Correlation between proton (antiproton) daughter MC pt and Xi/Omega MC pt (to apply Geant/Fluka correction)
      f2dHistAsMCptProtonMCptXiMinus(0),
      f2dHistAsMCptAntiprotonMCptXiPlus(0),
      f2dHistAsMCptProtonMCptOmegaMinus(0),
      f2dHistAsMCptAntiprotonMCptOmegaPlus(0),
      // -- Containers                       
      fCFContCascadePIDAsXiMinus(0),
      fCFContCascadePIDAsXiPlus(0),
      fCFContCascadePIDAsOmegaMinus(0),
      fCFContCascadePIDAsOmegaPlus(0),
      fCFContAsCascadeCuts(0)


    //____Costructor____
    {
      // Define input and output slots here
      // Input slot #0 works with a TChain
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
              
      DefineOutput(1, TList::Class());
      DefineOutput(2, AliCFContainer::Class());
      DefineOutput(3, AliCFContainer::Class());
      DefineOutput(4, AliCFContainer::Class());
      DefineOutput(5, AliCFContainer::Class());
      DefineOutput(6, AliCFContainer::Class());
    }

    //____Destructor____
    AliAnalysisTaskCheckPerformanceCascadepp::~AliAnalysisTaskCheckPerformanceCascadepp()
    {
      // For all TH1, 2, 3 HnSparse and CFContainer are in the fListCascade TList.
      // They will be deleted when fListCascade is deleted by the TSelector dtor
      // Because of TList::SetOwner()
      if (fListHistCascade && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())              {delete fListHistCascade;              fListHistCascade = 0x0;}  
      if (fCFContCascadePIDAsXiMinus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())    {delete fCFContCascadePIDAsXiMinus;    fCFContCascadePIDAsXiMinus = 0x0;}
      if (fCFContCascadePIDAsXiPlus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())     {delete fCFContCascadePIDAsXiPlus;     fCFContCascadePIDAsXiPlus = 0x0;}
      if (fCFContCascadePIDAsOmegaMinus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {delete fCFContCascadePIDAsOmegaMinus; fCFContCascadePIDAsOmegaMinus = 0x0;}
      if (fCFContCascadePIDAsOmegaPlus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())  {delete fCFContCascadePIDAsOmegaPlus;  fCFContCascadePIDAsOmegaPlus = 0x0;}
      if (fCFContAsCascadeCuts && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())          {delete fCFContAsCascadeCuts;          fCFContAsCascadeCuts = 0x0;}
      if (fESDtrackCuts)                                                                             {delete fESDtrackCuts;                 fESDtrackCuts = 0x0;}
    }


//________________________________________________________________________
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// - CREATE THE OUTPUTS
void AliAnalysisTaskCheckPerformanceCascadepp::UserCreateOutputObjects() {
  // Create histograms
  // Called once

 // - Option for AliLog: to suppress the extensive info prompted by a run with MC
 AliLog::SetGlobalLogLevel(AliLog::kError); 

 // - Definition of the output datamembers	
 fListHistCascade = new TList();
 fListHistCascade->SetOwner(); // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

 //-----------------------------------------------
 // Particle Identification Setup (new PID object)
 //-----------------------------------------------
 AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
 AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
 fPIDResponse = inputHandler->GetPIDResponse();
        
 // - Only used to get the number of primary reconstructed tracks
 if (! fESDtrackCuts ){
      fESDtrackCuts = new AliESDtrackCuts();
 }
 if (! fUtils){
      fUtils = new AliAnalysisUtils();
 }
 
 //----------------------
 // Initialize the histos
 //----------------------

 //----------------------------------
 // - General binning definition
 Double_t ptBinLimits[101];
 for (Int_t iptbin = 0; iptbin<101; ++iptbin) ptBinLimits[iptbin]=iptbin*0.1;
 Double_t yBinLimits[111];
 for (Int_t iybin = 0; iybin<111; ++iybin) yBinLimits[iybin]=-1.1+iybin*0.02;
 Double_t ctauBinLimits[112];
 for (Int_t ict = 0; ict<112; ++ict) ctauBinLimits[ict] = (Double_t) (ict-1.); 
 
 //------------------
 // - General plots
   // - Cascades multiplicity plots 
   if (! fHistCascadeMultiplicityBeforeAnySel) {
        fHistCascadeMultiplicityBeforeAnySel = new TH1F("fHistCascadeMultiplicityBeforeAnySel", "Cascades per event (before any selections);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityBeforeAnySel);
   }
   if (! fHistCascadeMultiplicityAfterSDDstatusSel) {
        fHistCascadeMultiplicityAfterSDDstatusSel = new TH1F("fHistCascadeMultiplicityAfterSDDstatusSel", "Cascades per event (after the SDD status selection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterSDDstatusSel);
   }
   if (! fHistCascadeMultiplicityAfterDAQincompleteEvRej) {
        fHistCascadeMultiplicityAfterDAQincompleteEvRej = new TH1F("fHistCascadeMultiplicityAfterDAQincompleteEvRej", "Cascades per event (after DAQ incomplete event rejection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterDAQincompleteEvRej);
   }
   if (! fHistCascadeMultiplicityAfterSPDclustervstrackletSel) {
        fHistCascadeMultiplicityAfterSPDclustervstrackletSel = new TH1F("fHistCascadeMultiplicityAfterSPDclustervstrackletSel", "Cascades per event (after background rejection based on SPD cluster vs tracklet);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterSPDclustervstrackletSel);
   }
   if (! fHistCascadeMultiplicityAfterPileupRej) {
        fHistCascadeMultiplicityAfterPileupRej = new TH1F("fHistCascadeMultiplicityAfterPileupRej", "Cascades per event (after pile-up events rejection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterPileupRej);
   }
   if (! fHistCascadeMultiplicityAfterPhysicsSel) {
        fHistCascadeMultiplicityAfterPhysicsSel = new TH1F("fHistCascadeMultiplicityAfterPhysicsSel", "Cascades per event (after physics selection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterPhysicsSel);
   }
   if (! fHistCascadeMultiplicityAfterRevertexing) {
        fHistCascadeMultiplicityAfterRevertexing = new TH1F("fHistCascadeMultiplicityAfterRevertexing", "Cascades per event (after physics selection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterRevertexing);
   }
   if (! fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel) {
        fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel = new TH1F("fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel", "Cascades per event (for selected events with well-established PV);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel);
   }
   if (! fHistCascadeMultiplicityAfterSPDresolution) {
        fHistCascadeMultiplicityAfterSPDresolution = new TH1F("fHistCascadeMultiplicityAfterSPDresolution", "Cascades per event (after physics selection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterSPDresolution);
   }
   if (! fHistCascadeMultiplicityAfterVerticesProximity) {
        fHistCascadeMultiplicityAfterVerticesProximity = new TH1F("fHistCascadeMultiplicityAfterVerticesProximity", "Cascades per event (after physics selection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterVerticesProximity);
   }
   if (! fHistCascadeMultiplicityAfterZprimVtxPosSel) {
        fHistCascadeMultiplicityAfterZprimVtxPosSel = new TH1F("fHistCascadeMultiplicityAfterZprimVtxPosSel", "Cascades per event (after vertex cut selection);Nbr of Cascades/Evt;Events", 20, 0, 20);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterZprimVtxPosSel);
   }
   // - Vertex position plots
   // -- After all the event selections but the |Zvtx| < 10 cm cut
   if (! fHistPVx ){
        fHistPVx = new TH1F("fHistPVx", "Best PV position in x; x (cm); Events", 2000, -0.5, 0.5);
        fListHistCascade->Add(fHistPVx);
   }
   if (! fHistPVy ){
        fHistPVy = new TH1F("fHistPVy", "Best PV position in y; y (cm); Events", 4000, -1.0, 1.0);
        fListHistCascade->Add(fHistPVy);
   }
   if (! fHistPVz ){
        fHistPVz = new TH1F("fHistPVz", "Best PV position in z; z (cm); Events", 400, -20, 20);
        fListHistCascade->Add(fHistPVz);
   }
   // -- After all the event selections
   if (! fHistPVxAnalysis ){
        fHistPVxAnalysis = new TH1F("fHistPVxAnalysis", "Best PV position in x (after events selections); x (cm); Events", 2000, -0.5, 0.5);
        fListHistCascade->Add(fHistPVxAnalysis);
   }
   if (! fHistPVyAnalysis ){
        fHistPVyAnalysis = new TH1F("fHistPVyAnalysis", "Best PV position in y (after events selections); y (cm); Events" , 4000, -1.0, 1.0);
        fListHistCascade->Add(fHistPVyAnalysis);
   }
   if (! fHistPVzAnalysis ){
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

 //--------------------------
 // - Generated cascade plots
   // - Generated cascade multiplicity distributions (for single event)
   if (! fHistnXiPlusPerEvTot){
        fHistnXiPlusPerEvTot = new TH1F("fHistnXiPlusPerEvTot", "", 25, 0, 25);
        fListHistCascade->Add(fHistnXiPlusPerEvTot);
   }
   if (! fHistnXiMinusPerEvTot){
        fHistnXiMinusPerEvTot = new TH1F("fHistnXiMinusPerEvTot", "", 25, 0, 25);
        fListHistCascade->Add(fHistnXiMinusPerEvTot);
   }
   if (! fHistnOmegaPlusPerEvTot){
        fHistnOmegaPlusPerEvTot = new TH1F("fHistnOmegaPlusPerEvTot", "", 25, 0, 25);
        fListHistCascade->Add(fHistnOmegaPlusPerEvTot);
   }
   if (! fHistnOmegaMinusPerEvTot){
        fHistnOmegaMinusPerEvTot = new TH1F("fHistnOmegaMinusPerEvTot", "", 25, 0, 25);
        fListHistCascade->Add(fHistnOmegaMinusPerEvTot);   
   }
   if (! fHistnXiPlusPerEv){
        fHistnXiPlusPerEv = new TH1F("fHistnXiPlusPerEv", "", 25, 0, 25);
        fListHistCascade->Add(fHistnXiPlusPerEv);
   }
   if (! fHistnXiMinusPerEv){
        fHistnXiMinusPerEv = new TH1F("fHistnXiMinusPerEv", "", 25, 0, 25);
        fListHistCascade->Add(fHistnXiMinusPerEv);
   }
   if (! fHistnOmegaPlusPerEv){
        fHistnOmegaPlusPerEv = new TH1F("fHistnOmegaPlusPerEv", "", 25, 0, 25);
        fListHistCascade->Add(fHistnOmegaPlusPerEv);
   }
   if (! fHistnOmegaMinusPerEv){
        fHistnOmegaMinusPerEv = new TH1F("fHistnOmegaMinusPerEv", "", 25, 0, 25);
        fListHistCascade->Add(fHistnOmegaMinusPerEv);
   }
   // - Efficiency denominator plots
   // -- (A) Before any event selection
   if (! f3dHistGenPtVsGenYvsNtracksXiMinus_A) {
        f3dHistGenPtVsGenYvsNtracksXiMinus_A = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinus_A", "MC P_{t} Vs MC Y of Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinus_A);
   }
   if (! f3dHistGenPtVsGenYvsNtracksXiPlus_A) {
        f3dHistGenPtVsGenYvsNtracksXiPlus_A = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlus_A", "MC P_{t} Vs MC Y of Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlus_A);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaMinus_A) {
        f3dHistGenPtVsGenYvsNtracksOmegaMinus_A = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinus_A", "MC P_{t} Vs MC Y of Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinus_A);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaPlus_A) {
        f3dHistGenPtVsGenYvsNtracksOmegaPlus_A = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlus_A", "MC P_{t} Vs MC Y of Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlus_A);
   }
   // -- (B) After preliminary event selections (IsIncompleteDAQ, Tracklet vs Clusters Cut, SPD Pileup 
   if (! f3dHistGenPtVsGenYvsNtracksXiMinus_B) {
        f3dHistGenPtVsGenYvsNtracksXiMinus_B = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinus_B", "MC P_{t} Vs MC Y of Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinus_B);
   }
   if (! f3dHistGenPtVsGenYvsNtracksXiPlus_B) {
        f3dHistGenPtVsGenYvsNtracksXiPlus_B = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlus_B", "MC P_{t} Vs MC Y of Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlus_B);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaMinus_B) {
        f3dHistGenPtVsGenYvsNtracksOmegaMinus_B = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinus_B", "MC P_{t} Vs MC Y of Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinus_B);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaPlus_B) {
        f3dHistGenPtVsGenYvsNtracksOmegaPlus_B = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlus_B", "MC P_{t} Vs MC Y of Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlus_B);
   }
   if (! f3dHistGenPtVsGenctauvsYXiMinus_B) {
        f3dHistGenPtVsGenctauvsYXiMinus_B = new TH3D("f3dHistGenPtVsGenctauvsYXiMinus_B", "MC P_{t} Vs MC ctau Vs Y of Gen #Xi^{-}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
        fListHistCascade->Add(f3dHistGenPtVsGenctauvsYXiMinus_B);
   }
   if (! f3dHistGenPtVsGenctauvsYXiPlus_B) {
        f3dHistGenPtVsGenctauvsYXiPlus_B = new TH3D("f3dHistGenPtVsGenctauvsYXiPlus_B", "MC P_{t} Vs MC ctau Vs Y of Gen #Xi^{+}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
        fListHistCascade->Add(f3dHistGenPtVsGenctauvsYXiPlus_B);
   }
   if (! f3dHistGenPtVsGenctauvsYOmegaMinus_B) {
        f3dHistGenPtVsGenctauvsYOmegaMinus_B = new TH3D("f3dHistGenPtVsGenctauvsYOmegaMinus_B", "MC P_{t} Vs MC ctau Vs Y of Gen #Omega^{-}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
        fListHistCascade->Add(f3dHistGenPtVsGenctauvsYOmegaMinus_B);
   }
   if (! f3dHistGenPtVsGenctauvsYOmegaPlus_B) {
        f3dHistGenPtVsGenctauvsYOmegaPlus_B = new TH3D("f3dHistGenPtVsGenctauvsYOmegaPlus_B", "MC P_{t} Vs MC ctau Vs Y of Gen #Omega^{+}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
        fListHistCascade->Add(f3dHistGenPtVsGenctauvsYOmegaPlus_B);
   }
   // -- (C) After physics selection
   if (! f3dHistGenPtVsGenYvsNtracksXiMinus_C) {
        f3dHistGenPtVsGenYvsNtracksXiMinus_C = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinus_C", "MC P_{t} Vs MC Y of Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinus_C);
   }
   if (! f3dHistGenPtVsGenYvsNtracksXiPlus_C) {
        f3dHistGenPtVsGenYvsNtracksXiPlus_C = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlus_C", "MC P_{t} Vs MC Y of Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlus_C);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaMinus_C) {
        f3dHistGenPtVsGenYvsNtracksOmegaMinus_C = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinus_C", "MC P_{t} Vs MC Y of Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinus_C);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaPlus_C) {
        f3dHistGenPtVsGenYvsNtracksOmegaPlus_C = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlus_C", "MC P_{t} Vs MC Y of Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlus_C);
   }
   // -- (D) After re-vertexing
   if (! f3dHistGenPtVsGenYvsNtracksXiMinus_D) {
        f3dHistGenPtVsGenYvsNtracksXiMinus_D = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinus_D", "MC P_{t} Vs MC Y of Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinus_D);
   }
   if (! f3dHistGenPtVsGenYvsNtracksXiPlus_D) {
        f3dHistGenPtVsGenYvsNtracksXiPlus_D = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlus_D", "MC P_{t} Vs MC Y of Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlus_D);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaMinus_D) {
        f3dHistGenPtVsGenYvsNtracksOmegaMinus_D = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinus_D", "MC P_{t} Vs MC Y of Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinus_D);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaPlus_D) {
        f3dHistGenPtVsGenYvsNtracksOmegaPlus_D = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlus_D", "MC P_{t} Vs MC Y of Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlus_D);
   }
   // -- (E)  
   if (! f3dHistGenPtVsGenYvsNtracksXiMinus_E) {
        f3dHistGenPtVsGenYvsNtracksXiMinus_E = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinus_E", "MC P_{t} Vs MC Y of Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinus_E);
   }
   if (! f3dHistGenPtVsGenYvsNtracksXiPlus_E) {
        f3dHistGenPtVsGenYvsNtracksXiPlus_E = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlus_E", "MC P_{t} Vs MC Y of Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlus_E);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaMinus_E) {
        f3dHistGenPtVsGenYvsNtracksOmegaMinus_E = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinus_E", "MC P_{t} Vs MC Y of Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinus_E);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaPlus_E) {
        f3dHistGenPtVsGenYvsNtracksOmegaPlus_E = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlus_E", "MC P_{t} Vs MC Y of Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlus_E);
   }
   // -- (F)  
   if (! f3dHistGenPtVsGenYvsNtracksXiMinus_F) {
        f3dHistGenPtVsGenYvsNtracksXiMinus_F = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinus_F", "MC P_{t} Vs MC Y of Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinus_F);
   }
   if (! f3dHistGenPtVsGenYvsNtracksXiPlus_F) {
        f3dHistGenPtVsGenYvsNtracksXiPlus_F = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlus_F", "MC P_{t} Vs MC Y of Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlus_F);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaMinus_F) {
        f3dHistGenPtVsGenYvsNtracksOmegaMinus_F = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinus_F", "MC P_{t} Vs MC Y of Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinus_F);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaPlus_F) {
        f3dHistGenPtVsGenYvsNtracksOmegaPlus_F = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlus_F", "MC P_{t} Vs MC Y of Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlus_F);
   }
   // -- (G)  
   if (! f3dHistGenPtVsGenYvsNtracksXiMinus_G) {
        f3dHistGenPtVsGenYvsNtracksXiMinus_G = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinus_G", "MC P_{t} Vs MC Y of Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinus_G);
   }
   if (! f3dHistGenPtVsGenYvsNtracksXiPlus_G) {
        f3dHistGenPtVsGenYvsNtracksXiPlus_G = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlus_G", "MC P_{t} Vs MC Y of Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlus_G);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaMinus_G) {
        f3dHistGenPtVsGenYvsNtracksOmegaMinus_G = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinus_G", "MC P_{t} Vs MC Y of Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinus_G);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaPlus_G) {
        f3dHistGenPtVsGenYvsNtracksOmegaPlus_G = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlus_G", "MC P_{t} Vs MC Y of Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlus_G);
   }
   // -- (H)
   if (! f3dHistGenPtVsGenYvsNtracksXiMinus_H) {
        f3dHistGenPtVsGenYvsNtracksXiMinus_H = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinus_H", "MC P_{t} Vs MC Y of Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinus_H);
   }
   if (! f3dHistGenPtVsGenYvsNtracksXiPlus_H) {
        f3dHistGenPtVsGenYvsNtracksXiPlus_H = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlus_H", "MC P_{t} Vs MC Y of Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlus_H);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaMinus_H) {
        f3dHistGenPtVsGenYvsNtracksOmegaMinus_H = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinus_H", "MC P_{t} Vs MC Y of Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinus_H);
   }
   if (! f3dHistGenPtVsGenYvsNtracksOmegaPlus_H) {
        f3dHistGenPtVsGenYvsNtracksOmegaPlus_H = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlus_H", "MC P_{t} Vs MC Y of Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
        fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlus_H);
   }
   if (! f3dHistGenPtVsGenctauvsYXiMinus_H) {
        f3dHistGenPtVsGenctauvsYXiMinus_H = new TH3D("f3dHistGenPtVsGenctauvsYXiMinus_H", "MC P_{t} Vs MC ctau Vs Y of Gen #Xi^{-}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
        fListHistCascade->Add(f3dHistGenPtVsGenctauvsYXiMinus_H);
   }
   if (! f3dHistGenPtVsGenctauvsYXiPlus_H) {
        f3dHistGenPtVsGenctauvsYXiPlus_H = new TH3D("f3dHistGenPtVsGenctauvsYXiPlus_H", "MC P_{t} Vs MC ctau Vs Y of Gen #Xi^{+}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
        fListHistCascade->Add(f3dHistGenPtVsGenctauvsYXiPlus_H);
   }
   if (! f3dHistGenPtVsGenctauvsYOmegaMinus_H) {
        f3dHistGenPtVsGenctauvsYOmegaMinus_H = new TH3D("f3dHistGenPtVsGenctauvsYOmegaMinus_H", "MC P_{t} Vs MC ctau Vs Y of Gen #Omega^{-}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
        fListHistCascade->Add(f3dHistGenPtVsGenctauvsYOmegaMinus_H);
   }
   if (! f3dHistGenPtVsGenctauvsYOmegaPlus_H) {
        f3dHistGenPtVsGenctauvsYOmegaPlus_H = new TH3D("f3dHistGenPtVsGenctauvsYOmegaPlus_H", "MC P_{t} Vs MC ctau Vs Y of Gen #Omega^{+}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
        fListHistCascade->Add(f3dHistGenPtVsGenctauvsYOmegaPlus_H);
   }
   // - Many observable distributions for mother and daugthers: pseudo-rapidity, theta, Pt vs Y, Pt
   // -- xi minus
   if (! fHistEtaGenCascXiMinus) {
        fHistEtaGenCascXiMinus = new TH1F("fHistEtaGenCascXiMinus", "#eta of any gen. #Xi^{-}; #eta; Number of Casc", 200, -10, 10);
        fListHistCascade->Add(fHistEtaGenCascXiMinus);
   }
   if (! fHistThetaGenCascXiMinus) {
        fHistThetaGenCascXiMinus = new TH1F("fHistThetaGenCascXiMinus", "#theta of gen. #Xi^{-}; #theta; Number of Casc.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaGenCascXiMinus);
   }
   if (! f2dHistGenPtVsGenYFdblXiMinus) {
        f2dHistGenPtVsGenYFdblXiMinus = new TH2D("f2dHistGenPtVsGenYFdblXiMinus", "MC P_{t} Vs MC Y of findable Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistGenPtVsGenYFdblXiMinus);
   }
   if (! fHistThetaLambdaXiMinus) {
        fHistThetaLambdaXiMinus = new TH1F("fHistThetaLambdaXiMinus", "#theta of gen. #Lambda (Xi dghter); #theta_{#Lambda}; Number of #Lambda^0", 200, -10, 190);
        fListHistCascade->Add(fHistThetaLambdaXiMinus);
   }
   if (! fHistThetaBachXiMinus) {
        fHistThetaBachXiMinus = new TH1F("fHistThetaBachXiMinus", "#theta of gen. Bach.; #theta_{Bach}; Number of Bach.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaBachXiMinus);
   }
   if (! fHistThetaMesDghterXiMinus) {
        fHistThetaMesDghterXiMinus = new TH1F("fHistThetaMesDghterXiMinus", "#theta of gen. Meson #Lambda dghter; #theta_{MesDght}; Number of Mes.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaMesDghterXiMinus);
   }
   if (! fHistThetaBarDghterXiMinus) {
        fHistThetaBarDghterXiMinus = new TH1F("fHistThetaBarDghterXiMinus", "#theta of gen. Baryon #Lambda dghter; #theta_{BarDght}; Number of Bar.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaBarDghterXiMinus);
   }
   if (! fHistEtaLambdaXiMinus) {
        fHistEtaLambdaXiMinus = new TH1F("fHistEtaLambdaXiMinus", "#eta of gen. #Lambda (Xi dghter); #eta_{#Lambda}; Number of #Lambda^0", 200, -10, 10);
        fListHistCascade->Add(fHistEtaLambdaXiMinus);
   }
   if (! fHistEtaBachXiMinus) {
        fHistEtaBachXiMinus = new TH1F("fHistEtaBachXiMinus", "#eta of gen. Bach.; #eta_{Bach}; Number of Bach.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaBachXiMinus);
   }
   if (! fHistEtaMesDghterXiMinus) {
        fHistEtaMesDghterXiMinus = new TH1F("fHistEtaMesDghterXiMinus", "#eta of gen. Meson #Lambda dghter; #eta_{MesDght}; Number of Mes.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaMesDghterXiMinus);
   }
   if (! fHistEtaBarDghterXiMinus) {
        fHistEtaBarDghterXiMinus = new TH1F("fHistEtaBarDghterXiMinus", "#eta of gen. Baryon #Lambda dghter; #eta_{BarDght}; Number of Bar.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaBarDghterXiMinus);
   }
   if (! fHistPtBachXiMinus) {
        fHistPtBachXiMinus = new TH1F("fHistPtBachXiMinus", "p_{t} of gen. Bach.; pt_{Bach}; Number of Bach.", 200, 0, 10);
        fListHistCascade->Add(fHistPtBachXiMinus);
   }
   if (! fHistPtMesDghterXiMinus) {
        fHistPtMesDghterXiMinus = new TH1F("fHistPtMesDghterXiMinus", "p_{t} of gen. Meson #Lambda dghter; pt_{MesDght}; Number of Mes.", 200, 0, 10);
        fListHistCascade->Add(fHistPtMesDghterXiMinus);
   }
   if (! fHistPtBarDghterXiMinus) {
        fHistPtBarDghterXiMinus = new TH1F("fHistPtBarDghterXiMinus", "p_{t} of gen. Baryon #Lambda dghter; pt_{BarDght}; Number of Bar.", 200, 0, 10);
        fListHistCascade->Add(fHistPtBarDghterXiMinus);
   }
   // -- xi plus
   if (! fHistEtaGenCascXiPlus) {
        fHistEtaGenCascXiPlus = new TH1F("fHistEtaGenCascXiPlus", "#eta of any gen. #Xi^{+}; #eta; Number of Casc", 200, -10, 10);
        fListHistCascade->Add(fHistEtaGenCascXiPlus);
   }
   if (! fHistThetaGenCascXiPlus) {
        fHistThetaGenCascXiPlus = new TH1F("fHistThetaGenCascXiPlus", "#theta of gen. #Xi^{+}; #theta; Number of Casc.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaGenCascXiPlus);
   }
   if (! f2dHistGenPtVsGenYFdblXiPlus) {
        f2dHistGenPtVsGenYFdblXiPlus = new TH2D("f2dHistGenPtVsGenYFdblXiPlus", "MC P_{t} Vs MC Y of findable Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistGenPtVsGenYFdblXiPlus);
   }
   if (! fHistThetaLambdaXiPlus) {
        fHistThetaLambdaXiPlus = new TH1F("fHistThetaLambdaXiPlus", "#theta of gen. #Lambda (Xi dghter); #theta_{#Lambda}; Number of #Lambda", 200, -10, 190);
        fListHistCascade->Add(fHistThetaLambdaXiPlus);
   }
   if (! fHistThetaBachXiPlus) {
        fHistThetaBachXiPlus = new TH1F("fHistThetaBachXiPlus", "#theta of gen. Bach.; #theta_{Bach}; Number of Bach.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaBachXiPlus);
   }
   if (! fHistThetaMesDghterXiPlus) {
        fHistThetaMesDghterXiPlus = new TH1F("fHistThetaMesDghterXiPlus", "#theta of gen. Meson #Lambda dghter; #theta_{MesDght}; Number of Mes.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaMesDghterXiPlus);
   }
   if (! fHistThetaBarDghterXiPlus) {
        fHistThetaBarDghterXiPlus = new TH1F("fHistThetaBarDghterXiPlus", "#theta of gen. Baryon #Lambda dghter; #theta_{BarDght}; Number of Bar.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaBarDghterXiPlus);
   }
   if (! fHistEtaLambdaXiPlus) {
        fHistEtaLambdaXiPlus = new TH1F("fHistEtaLambdaXiPlus", "#eta of gen. #Lambda (Xi dghter); #eta_{#Lambda}; Number of #Lambda", 200, -10, 10);
        fListHistCascade->Add(fHistEtaLambdaXiPlus);
   }
   if (! fHistEtaBachXiPlus) {
        fHistEtaBachXiPlus = new TH1F("fHistEtaBachXiPlus", "#eta of gen. Bach.; #eta_{Bach}; Number of Bach.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaBachXiPlus);
   }
   if (! fHistEtaMesDghterXiPlus) {
        fHistEtaMesDghterXiPlus = new TH1F("fHistEtaMesDghterXiPlus", "#eta of gen. Meson #Lambda dghter; #eta_{MesDght}; Number of Mes.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaMesDghterXiPlus);
   }
   if (! fHistEtaBarDghterXiPlus) {
        fHistEtaBarDghterXiPlus = new TH1F("fHistEtaBarDghterXiPlus", "#eta of gen. Baryon #Lambda dghter; #eta_{BarDght}; Number of Bar.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaBarDghterXiPlus);
   }
   if (! fHistPtBachXiPlus) {
        fHistPtBachXiPlus = new TH1F("fHistPtBachXiPlus", "p_{t} of gen. Bach.; pt_{Bach}; Number of Bach.", 200, 0, 10);
        fListHistCascade->Add(fHistPtBachXiPlus);
   }
   if (! fHistPtMesDghterXiPlus) {
        fHistPtMesDghterXiPlus = new TH1F("fHistPtMesDghterXiPlus", "p_{t} of gen. Meson #Lambda dghter; pt_{MesDght}; Number of Mes.", 200, 0, 10);
        fListHistCascade->Add(fHistPtMesDghterXiPlus);
   }
   if (! fHistPtBarDghterXiPlus) {
        fHistPtBarDghterXiPlus = new TH1F("fHistPtBarDghterXiPlus", "p_{t} of gen. Baryon #Lambda dghter); pt_{BarDght}; Number of Bar.", 200, 0, 10);
        fListHistCascade->Add(fHistPtBarDghterXiPlus);
   }
   // -- omega minus
   if (! fHistEtaGenCascOmegaMinus) {
        fHistEtaGenCascOmegaMinus = new TH1F("fHistEtaGenCascOmegaMinus", "#eta of any gen. #Omega^{-}; #eta; Number of Casc", 200, -10, 10);
        fListHistCascade->Add(fHistEtaGenCascOmegaMinus);
   }
   if (! fHistThetaGenCascOmegaMinus) {
        fHistThetaGenCascOmegaMinus = new TH1F("fHistThetaGenCascOmegaMinus", "#theta of gen. #Omega^{-}; #theta; Number of Casc.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaGenCascOmegaMinus);
   }
   if (! f2dHistGenPtVsGenYFdblOmegaMinus) {
        f2dHistGenPtVsGenYFdblOmegaMinus = new TH2D("f2dHistGenPtVsGenYFdblOmegaMinus", "MC P_{t} Vs MC Y of findable Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistGenPtVsGenYFdblOmegaMinus);
   }
   if (! fHistThetaLambdaOmegaMinus) {
        fHistThetaLambdaOmegaMinus = new TH1F("fHistThetaLambdaOmegaMinus", "#theta of gen. #Lambda (Omega dghter); #theta_{#Lambda}; Number of #Lambda", 200, -10, 190);
        fListHistCascade->Add(fHistThetaLambdaOmegaMinus);
   }
   if (! fHistThetaBachOmegaMinus) {
        fHistThetaBachOmegaMinus = new TH1F("fHistThetaBachOmegaMinus", "#theta of gen. Bach.;#theta_{Bach};Number of Bach.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaBachOmegaMinus);
   }
   if (! fHistThetaMesDghterOmegaMinus) {
        fHistThetaMesDghterOmegaMinus = new TH1F("fHistThetaMesDghterOmegaMinus", "#theta of gen. Meson #Lambda dghter; #theta_{MesDght}; Number of Mes.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaMesDghterOmegaMinus);
   }
   if (! fHistThetaBarDghterOmegaMinus) {
        fHistThetaBarDghterOmegaMinus = new TH1F("fHistThetaBarDghterOmegaMinus", "#theta of gen. Baryon #Lambda dghter; #theta_{BarDght}; Number of Bar.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaBarDghterOmegaMinus);
   }
   if (! fHistEtaLambdaOmegaMinus) {
        fHistEtaLambdaOmegaMinus = new TH1F("fHistEtaLambdaOmegaMinus", "#eta of gen. #Lambda (Omega dghter); #eta_{#Lambda}; Number of #Lambda", 200, -10, 10);
        fListHistCascade->Add(fHistEtaLambdaOmegaMinus);
   }
   if (! fHistEtaBachOmegaMinus) {
        fHistEtaBachOmegaMinus = new TH1F("fHistEtaBachOmegaMinus", "#eta of gen. Bach.;#eta_{Bach};Number of Bach.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaBachOmegaMinus);
   }
   if (! fHistEtaMesDghterOmegaMinus) {
        fHistEtaMesDghterOmegaMinus = new TH1F("fHistEtaMesDghterOmegaMinus", "#eta of gen. Meson #Lambda dghter; #eta_{MesDght}; Number of Mes.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaMesDghterOmegaMinus);
   }
   if (! fHistEtaBarDghterOmegaMinus) {
        fHistEtaBarDghterOmegaMinus = new TH1F("fHistEtaBarDghterOmegaMinus", "#eta of gen. Baryon #Lambda dghter; #eta_{BarDght}; Number of Bar.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaBarDghterOmegaMinus);
   }
   if (! fHistPtBachOmegaMinus) {
        fHistPtBachOmegaMinus = new TH1F("fHistPtBachOmegaMinus", "p_{t} of gen. Bach.; pt_{Bach}; Number of Bach.", 200, 0, 10);
        fListHistCascade->Add(fHistPtBachOmegaMinus);
   }
   if (! fHistPtMesDghterOmegaMinus) {
        fHistPtMesDghterOmegaMinus = new TH1F("fHistPtMesDghterOmegaMinus", "p_{t} of gen. Meson #Lambda dghter); pt_{MesDght}; Number of Mes.", 200, 0, 10);
        fListHistCascade->Add(fHistPtMesDghterOmegaMinus);
   }
   if (! fHistPtBarDghterOmegaMinus) {
        fHistPtBarDghterOmegaMinus = new TH1F("fHistPtBarDghterOmegaMinus", "p_{t} of gen. Baryon #Lambda dghter); pt_{BarDght}; Number of Bar.", 200, 0, 10);
        fListHistCascade->Add(fHistPtBarDghterOmegaMinus);
   }
   // -- omega plus
   if (! fHistEtaGenCascOmegaPlus) {
        fHistEtaGenCascOmegaPlus = new TH1F("fHistEtaGenCascOmegaPlus", "#eta of any gen. #Omega^{+}; #eta; Number of Casc", 200, -10, 10);
        fListHistCascade->Add(fHistEtaGenCascOmegaPlus);
   }
   if (! fHistThetaGenCascOmegaPlus) {
        fHistThetaGenCascOmegaPlus = new TH1F("fHistThetaGenCascOmegaPlus", "#theta of gen. #Omega^{+}; #theta; Number of Casc.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaGenCascOmegaPlus);
   }
   if (! f2dHistGenPtVsGenYFdblOmegaPlus) {
        f2dHistGenPtVsGenYFdblOmegaPlus = new TH2D("f2dHistGenPtVsGenYFdblOmegaPlus", "MC P_{t} Vs MC Y of findable Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistGenPtVsGenYFdblOmegaPlus);
   }
   if (! fHistThetaLambdaOmegaPlus) {
        fHistThetaLambdaOmegaPlus = new TH1F("fHistThetaLambdaOmegaPlus", "#theta of gen. #Lambda (Omega dghter); #theta_{#Lambda}; Number of #Lambda", 200, -10, 190);
        fListHistCascade->Add(fHistThetaLambdaOmegaPlus);
   }
   if (! fHistThetaBachOmegaPlus) {
        fHistThetaBachOmegaPlus = new TH1F("fHistThetaBachOmegaPlus", "#theta of gen. Bach.; #theta_{Bach}; Number of Bach.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaBachOmegaPlus);
   }
   if (! fHistThetaMesDghterOmegaPlus) {
        fHistThetaMesDghterOmegaPlus = new TH1F("fHistThetaMesDghterOmegaPlus", "#theta of gen. Meson #Lambda dghter; #theta_{MesDght}; Number of Mes.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaMesDghterOmegaPlus);
   }
   if (! fHistThetaBarDghterOmegaPlus) {
        fHistThetaBarDghterOmegaPlus = new TH1F("fHistThetaBarDghterOmegaPlus", "#theta of gen. Baryon #Lambda dghter; #theta_{BarDght}; Number of Bar.", 200, -10, 190);
        fListHistCascade->Add(fHistThetaBarDghterOmegaPlus);
   }
   if (! fHistEtaLambdaOmegaPlus) {
        fHistEtaLambdaOmegaPlus = new TH1F("fHistEtaLambdaOmegaPlus", "#eta of gen. #Lambda (Omega dghter); #eta_{#Lambda}; Number of #Lambda", 200, -10, 10);
        fListHistCascade->Add(fHistEtaLambdaOmegaPlus);
   }
   if (! fHistEtaBachOmegaPlus) {
        fHistEtaBachOmegaPlus = new TH1F("fHistEtaBachOmegaPlus", "#eta of gen. Bach.; #eta_{Bach}; Number of Bach.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaBachOmegaPlus);
   }
   if (! fHistEtaMesDghterOmegaPlus) {
        fHistEtaMesDghterOmegaPlus = new TH1F("fHistEtaMesDghterOmegaPlus", "#eta of gen. Meson #Lambda dghter; #eta_{MesDght}; Number of Mes.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaMesDghterOmegaPlus);
   }
   if (! fHistEtaBarDghterOmegaPlus) {
        fHistEtaBarDghterOmegaPlus = new TH1F("fHistEtaBarDghterOmegaPlus", "#eta of gen. Baryon #Lambda dghter; #eta_{BarDght}; Number of Bar.", 200, -10, 10);
        fListHistCascade->Add(fHistEtaBarDghterOmegaPlus);
   }
   if (! fHistPtBachOmegaPlus) {
        fHistPtBachOmegaPlus = new TH1F("fHistPtBachOmegaPlus", "p_{t} of gen. Bach.; pt_{Bach}; Number of Bach.", 200, 0, 10);
        fListHistCascade->Add(fHistPtBachOmegaPlus);
   }
   if (! fHistPtMesDghterOmegaPlus) {
        fHistPtMesDghterOmegaPlus = new TH1F("fHistPtMesDghterOmegaPlus", "p_{t} of gen. Meson #Lambda dghter; pt_{MesDght}; Number of Mes.", 200, 0, 10);
        fListHistCascade->Add(fHistPtMesDghterOmegaPlus);
   }
   if (! fHistPtBarDghterOmegaPlus) {
        fHistPtBarDghterOmegaPlus = new TH1F("fHistPtBarDghterOmegaPlus", "p_{t} of gen. Baryon #Lambda dghter); pt_{BarDght}; Number of Bar.", 200, 0, 10);
        fListHistCascade->Add(fHistPtBarDghterOmegaPlus);
   }
 //-------------------------------------------------------------------------
 // - Any reconstructed cascades + reconstructed cascades associated with MC
   // - Associated cascade multiplicity distributions (for single event)
   if (! fHistnAssoXiMinus) {
        fHistnAssoXiMinus= new TH1F("fHistnAssoXiMinus", "", 25, 0, 25);
        fListHistCascade->Add(fHistnAssoXiMinus);
   }
   if (! fHistnAssoXiPlus) {
        fHistnAssoXiPlus= new TH1F("fHistnAssoXiPlus", "", 25, 0, 25);
        fListHistCascade->Add(fHistnAssoXiPlus);
   }
   if (! fHistnAssoOmegaMinus){
        fHistnAssoOmegaMinus= new TH1F("fHistnAssoOmegaMinus", "", 25, 0, 25);
        fListHistCascade->Add(fHistnAssoOmegaMinus);
   }
   if (! fHistnAssoOmegaPlus) {
        fHistnAssoOmegaPlus= new TH1F("fHistnAssoOmegaPlus", "", 25, 0, 25);
        fListHistCascade->Add(fHistnAssoOmegaPlus);
   }
   // - PID Probability versus MC Pt(bachelor track)
   if (! f2dHistPIDprobaKaonVsMCPtBach ){
        f2dHistPIDprobaKaonVsMCPtBach  = new TH2F("f2dHistPIDprobaKaonVsMCPtBach", "Comb. PID proba to be K^{#pm} Vs MC Bach. Pt; Pt_{MC}(Bach.) (GeV/c); Comb. PID Proba (Bach. = K^{#pm})", 100, 0.0, 5.0, 110, 0.0, 1.10);
        fListHistCascade->Add(f2dHistPIDprobaKaonVsMCPtBach);
   }
   if (! f2dHistPIDprobaPionVsMCPtBach ){
        f2dHistPIDprobaPionVsMCPtBach  = new TH2F("f2dHistPIDprobaPionVsMCPtBach", "Comb. PID proba to be #pi^{#pm} Vs MC Bach. Pt; Pt_{MC}(Bach.) (GeV/c); Comb. PID Proba (Bach. = #pi^{#pm})", 100, 0.0, 5.0, 110, 0.0, 1.10);
        fListHistCascade->Add(f2dHistPIDprobaPionVsMCPtBach);
   }
   // - Invariant mass distribution for cascades candidates ASSOCIATED with MC
   if (! fHistAsMCMassXiMinus) {
        fHistAsMCMassXiMinus = new TH1F("fHistAsMCMassXiMinus", "#Xi^{-} candidates associated to MC; M( #Lambda , #pi^{-} ) (GeV/c^{2}); Counts", 400, 1.2, 2.0);
        fListHistCascade->Add(fHistAsMCMassXiMinus);
   }
   if (! fHistAsMCMassXiPlus) {
        fHistAsMCMassXiPlus = new TH1F("fHistAsMCMassXiPlus", "#Xi^{+} candidates associated to MC; M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2}); Counts", 400, 1.2, 2.0);
        fListHistCascade->Add(fHistAsMCMassXiPlus);
   }
   if (! fHistAsMCMassOmegaMinus) {
        fHistAsMCMassOmegaMinus = new TH1F("fHistAsMCMassOmegaMinus", "#Omega^{-} candidates associated to MC; M( #Lambda , K^{-} ) (GeV/c^{2}); Counts", 500, 1.5, 2.5);
        fListHistCascade->Add(fHistAsMCMassOmegaMinus);
   }
   if (! fHistAsMCMassOmegaPlus) {
        fHistAsMCMassOmegaPlus = new TH1F("fHistAsMCMassOmegaPlus", "#Omega^{+} candidates associated to MC; M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2}); Counts", 500, 1.5, 2.5);
        fListHistCascade->Add(fHistAsMCMassOmegaPlus);
   }
   // -  Generated Pt Vs generated Y of the cascade candidates associated with MC + having the proper maximum probability of combined PID for the bachelor
   if (! f2dHistAsMCandCombPIDGenPtVsGenYXiMinus) {
        f2dHistAsMCandCombPIDGenPtVsGenYXiMinus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYXiMinus", "MC P_{t} Vs MC Y of #Xi^{-} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYXiMinus);
   }
   if (! f2dHistAsMCandCombPIDGenPtVsGenYXiPlus) {
        f2dHistAsMCandCombPIDGenPtVsGenYXiPlus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYXiPlus", "MC P_{t} Vs MC Y of #Xi^{+} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYXiPlus);
   } 
   if (! f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus) {
        f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus", "MC P_{t} Vs MC Y of #Omega^{-} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus);
   }
   if (! f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus) {
        f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus", "MC P_{t} Vs MC Y of #Omega^{+} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus);
   }
   // - Generated Pt Vs Generated Y, for the cascade candidates associated with MC
   if (! f2dHistAsMCGenPtVsGenYXiMinus) {
        f2dHistAsMCGenPtVsGenYXiMinus = new TH2F("f2dHistAsMCGenPtVsGenYXiMinus", "MC P_{t} Vs MC Y of gen. #Xi^{-} (associated); Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistAsMCGenPtVsGenYXiMinus );
   }
   if (! f2dHistAsMCGenPtVsGenYXiPlus) {
        f2dHistAsMCGenPtVsGenYXiPlus = new TH2F("f2dHistAsMCGenPtVsGenYXiPlus", "MC P_{t} Vs MC Y of gen. #Xi^{+} (associated); Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistAsMCGenPtVsGenYXiPlus );
   }
   if (! f2dHistAsMCGenPtVsGenYOmegaMinus) {
        f2dHistAsMCGenPtVsGenYOmegaMinus = new TH2F("f2dHistAsMCGenPtVsGenYOmegaMinus", "MC P_{t} Vs MC Y of gen. #Omega^{-} (associated); Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistAsMCGenPtVsGenYOmegaMinus );
   }
   if (! f2dHistAsMCGenPtVsGenYOmegaPlus) {
        f2dHistAsMCGenPtVsGenYOmegaPlus = new TH2F("f2dHistAsMCGenPtVsGenYOmegaPlus", "MC P_{t} Vs MC Y of gen. #Omega^{+} (associated); Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
        fListHistCascade->Add(f2dHistAsMCGenPtVsGenYOmegaPlus );
   } 
   // - Generated Eta of the the cascade candidates associated with MC
   if (! fHistAsMCGenEtaXiMinus) {
        fHistAsMCGenEtaXiMinus = new TH1F("fHistAsMCGenEtaXiMinus", "#eta of gen. #Xi^{-} (associated); #eta; Count", 100, -5, 5);
        fListHistCascade->Add( fHistAsMCGenEtaXiMinus );
   }
   if (! fHistAsMCGenEtaXiPlus) {
        fHistAsMCGenEtaXiPlus = new TH1F("fHistAsMCGenEtaXiPlus", "#eta of gen. #Xi^{+} (associated); #eta; Count", 100, -5, 5);
        fListHistCascade->Add( fHistAsMCGenEtaXiPlus );
   }
   if (! fHistAsMCGenEtaOmegaMinus) {
        fHistAsMCGenEtaOmegaMinus = new TH1F("fHistAsMCGenEtaOmegaMinus", "#eta of gen. #Omega^{-} (associated);#eta;Number of Casc", 100, -5, 5);
        fListHistCascade->Add( fHistAsMCGenEtaOmegaMinus );
   }
   if (! fHistAsMCGenEtaOmegaPlus) {
        fHistAsMCGenEtaOmegaPlus = new TH1F("fHistAsMCGenEtaOmegaPlus", "#eta of gen. #Omega^{+} (associated); #eta; Count", 100, -5, 5);
        fListHistCascade->Add( fHistAsMCGenEtaOmegaPlus );
   }
   // - Resolution in Pt as function of generated Pt
   if (! f2dHistAsMCResPtXiMinus) {
        f2dHistAsMCResPtXiMinus = new TH2F("f2dHistAsMCResPtXiMinus", "Resolution in Pt reconstruction for #Xi^{-}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
        fListHistCascade->Add(f2dHistAsMCResPtXiMinus);
   }
   if (! f2dHistAsMCResPtXiPlus) {
        f2dHistAsMCResPtXiPlus = new TH2F("f2dHistAsMCResPtXiPlus", "Resolution in Pt reconstruction for #Xi^{+}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
        fListHistCascade->Add(f2dHistAsMCResPtXiPlus);
   }
   if (! f2dHistAsMCResPtOmegaMinus) {
        f2dHistAsMCResPtOmegaMinus = new TH2F("f2dHistAsMCResPtOmegaMinus", "Resolution in Pt reconstruction for #Omega^{-}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
        fListHistCascade->Add(f2dHistAsMCResPtOmegaMinus);
   }
   if (! f2dHistAsMCResPtOmegaPlus) {
        f2dHistAsMCResPtOmegaPlus = new TH2F("f2dHistAsMCResPtOmegaPlus", "Resolution in Pt reconstruction for #Omega^{+}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
        fListHistCascade->Add(f2dHistAsMCResPtOmegaPlus);
   }
   // - Resolution in R(2D) as function of generated R
   if (! f2dHistAsMCResRXiMinus) {
        f2dHistAsMCResRXiMinus = new TH2F("f2dHistAsMCResRXiMinus", "Resolution in transv. position for #Xi^{-}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
        fListHistCascade->Add(f2dHistAsMCResRXiMinus);
   }
   if (! f2dHistAsMCResRXiPlus) {
        f2dHistAsMCResRXiPlus = new TH2F("f2dHistAsMCResRXiPlus", "Resolution in transv. position for #Xi^{+}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
        fListHistCascade->Add(f2dHistAsMCResRXiPlus);
   } 
   if (! f2dHistAsMCResROmegaMinus) {
        f2dHistAsMCResROmegaMinus = new TH2F("f2dHistAsMCResROmegaMinus", "Resolution in transv. position for #Omega^{-}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
        fListHistCascade->Add(f2dHistAsMCResROmegaMinus);
   }
   if (! f2dHistAsMCResROmegaPlus) {
        f2dHistAsMCResROmegaPlus = new TH2F("f2dHistAsMCResROmegaPlus", "Resolution in transv. position for #Omega^{+}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
        fListHistCascade->Add(f2dHistAsMCResROmegaPlus);
   }
   // - Resolution in phi as function of generated Pt 
   if (! f2dHistAsMCResPhiXiMinus) {
        f2dHistAsMCResPhiXiMinus = new TH2F("f2dHistAsMCResPhiXiMinus", "Resolution in #phi for #Xi^{-}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);
        fListHistCascade->Add(f2dHistAsMCResPhiXiMinus);
   }
   if (! f2dHistAsMCResPhiXiPlus) {
        f2dHistAsMCResPhiXiPlus = new TH2F("f2dHistAsMCResPhiXiPlus", "Resolution in #phi for #Xi^{+}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);
        fListHistCascade->Add(f2dHistAsMCResPhiXiPlus);
   }
   if (! f2dHistAsMCResPhiOmegaMinus) {
        f2dHistAsMCResPhiOmegaMinus = new TH2F("f2dHistAsMCResPhiOmegaMinus", "Resolution in #phi for #Omega^{-}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);  
        fListHistCascade->Add(f2dHistAsMCResPhiOmegaMinus);
   }
   if (! f2dHistAsMCResPhiOmegaPlus) {
        f2dHistAsMCResPhiOmegaPlus = new TH2F("f2dHistAsMCResPhiOmegaPlus", "Resolution in #phi for #Omega^{+}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);
        fListHistCascade->Add(f2dHistAsMCResPhiOmegaPlus);
   }
   //  - Correlation between proton (antiproton) daughter MC pt and Xi/Omega MC pt (to apply Geant/Fluka correction)
   if (! f2dHistAsMCptProtonMCptXiMinus) {
        f2dHistAsMCptProtonMCptXiMinus = new TH2F("f2dHistAsMCptProtonMCptXiMinus", "Proton MC pt vs Xi- MC pt", 100, 0., 10., 100, 0., 10.); 
        fListHistCascade->Add(f2dHistAsMCptProtonMCptXiMinus);
   }
   if (! f2dHistAsMCptAntiprotonMCptXiPlus) {
        f2dHistAsMCptAntiprotonMCptXiPlus = new TH2F("f2dHistAsMCptAntiprotonMCptXiPlus", "Antiproton MC pt vs Xi+ MC pt", 100, 0., 10., 100, 0., 10.);
        fListHistCascade->Add(f2dHistAsMCptAntiprotonMCptXiPlus);
   }
   if (! f2dHistAsMCptProtonMCptOmegaMinus) {
        f2dHistAsMCptProtonMCptOmegaMinus = new TH2F("f2dHistAsMCptProtonMCptOmegaMinus", "Proton MC pt vs Omega- MC pt", 100, 0., 10., 100, 0., 10.);
        fListHistCascade->Add(f2dHistAsMCptProtonMCptOmegaMinus);
   }
   if (! f2dHistAsMCptAntiprotonMCptOmegaPlus) {
        f2dHistAsMCptAntiprotonMCptOmegaPlus = new TH2F("f2dHistAsMCptAntiprotonMCptOmegaPlus", "Antiproton MC pt vs Omega+ MC pt", 100, 0., 10., 100, 0., 10.);
        fListHistCascade->Add(f2dHistAsMCptAntiprotonMCptOmegaPlus);
   }

  // - CFContainer
  // - Usefull string
  TString sddstatus = "";
  if      (fCollidingSystem == 0 && fApplyEvSelSDDstatus && fwithSDD)  sddstatus = "_wSDDon";
  else if (fCollidingSystem == 0 && fApplyEvSelSDDstatus && !fwithSDD) sddstatus = "_wSDDoff";
  TString cfcontname_cascpidximinus = Form("fCFContCascadePIDAsXiMinus_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fMinTPCcrossrawoverfindable,fVtxRangeMax,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks);
  cfcontname_cascpidximinus.Append(Form("%s",sddstatus.Data()));
  cfcontname_cascpidximinus.Append(Form("%s",fSuffix.Data()));
  TString cfcontname_cascpidxiplus = Form("fCFContCascadePIDAsXiPlus_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fMinTPCcrossrawoverfindable,fVtxRangeMax,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks);
  cfcontname_cascpidxiplus.Append(Form("%s",sddstatus.Data()));
  cfcontname_cascpidxiplus.Append(Form("%s",fSuffix.Data()));
  TString cfcontname_cascpidomegaminus = Form("fCFContCascadePIDAsOmegaMinus_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fMinTPCcrossrawoverfindable,fVtxRangeMax,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks);
  cfcontname_cascpidomegaminus.Append(Form("%s",sddstatus.Data()));
  cfcontname_cascpidomegaminus.Append(Form("%s",fSuffix.Data()));
  TString cfcontname_cascpidomegaplus = Form("fCFContCascadePIDAsOmegaPlus_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fMinTPCcrossrawoverfindable,fVtxRangeMax,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks);
  cfcontname_cascpidomegaplus.Append(Form("%s",sddstatus.Data()));
  cfcontname_cascpidomegaplus.Append(Form("%s",fSuffix.Data()));
  TString cfcontname_casccuts = Form("fCFContAsCascadeCuts_minnTPCcls%i_clsfindratio%.1f_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fMinTPCcrossrawoverfindable,fVtxRangeMax,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks);
  cfcontname_casccuts.Append(Form("%s",sddstatus.Data()));
  cfcontname_casccuts.Append(Form("%s",fSuffix.Data()));
  // -- PID container Xi-
  if(! fCFContCascadePIDAsXiMinus)  {
     const Int_t  lNbSteps      =  7;
     const Int_t  lNbVariables  =  3;
       //Array for the number of bins in each dimension:
     Int_t lNbBinsPerVar[3] = {0};
     lNbBinsPerVar[0] = 100;
     lNbBinsPerVar[1] = 800;
     lNbBinsPerVar[2] = 22;
     fCFContCascadePIDAsXiMinus = new AliCFContainer(cfcontname_cascpidximinus,"Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
       //Setting the bin limits 
     fCFContCascadePIDAsXiMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
     fCFContCascadePIDAsXiMinus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
     fCFContCascadePIDAsXiMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
       //Setting the step title : one per PID case
     fCFContCascadePIDAsXiMinus->SetStepTitle(0, "No PID");
     fCFContCascadePIDAsXiMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
     fCFContCascadePIDAsXiMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
     fCFContCascadePIDAsXiMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
     fCFContCascadePIDAsXiMinus->SetStepTitle(4, "Comb. PID / Bachelor");
     fCFContCascadePIDAsXiMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
     fCFContCascadePIDAsXiMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
       //Setting the variable title, per axis
     fCFContCascadePIDAsXiMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
     fCFContCascadePIDAsXiMinus->SetVarTitle(1, "M( #Lambda , #pi^{-} ) (GeV/c^{2})");
     fCFContCascadePIDAsXiMinus->SetVarTitle(2, "Y_{#Xi}");
  }
  // -- PID container Xi+
  if(! fCFContCascadePIDAsXiPlus)  {
     const Int_t  lNbSteps      =  7;
     const Int_t  lNbVariables  =  3;
       //Array for the number of bins in each dimension :
     Int_t lNbBinsPerVar[3] = {0};
     lNbBinsPerVar[0] = 100;
     lNbBinsPerVar[1] = 800;
     lNbBinsPerVar[2] = 22;
     fCFContCascadePIDAsXiPlus = new AliCFContainer(cfcontname_cascpidxiplus,"Pt_{cascade} Vs M_{#bar{#Xi}^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
       //Setting the bin limits (valid  for v4-18-10-AN)
     fCFContCascadePIDAsXiPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
     fCFContCascadePIDAsXiPlus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
     fCFContCascadePIDAsXiPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity 
       //Setting the step title : one per PID case
     fCFContCascadePIDAsXiPlus->SetStepTitle(0, "No PID");
     fCFContCascadePIDAsXiPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
     fCFContCascadePIDAsXiPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
     fCFContCascadePIDAsXiPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
     fCFContCascadePIDAsXiPlus->SetStepTitle(4, "Comb. PID / Bachelor");
     fCFContCascadePIDAsXiPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
     fCFContCascadePIDAsXiPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");  
       //Setting the variable title, per axis
     fCFContCascadePIDAsXiPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
     fCFContCascadePIDAsXiPlus->SetVarTitle(1, "M( #Lambda , #pi^{+} ) (GeV/c^{2})");
     fCFContCascadePIDAsXiPlus->SetVarTitle(2, "Y_{#Xi}");
  }
  // -- PID container Omega-
  if(! fCFContCascadePIDAsOmegaMinus)  {
     const Int_t  lNbSteps      =  7;
     const Int_t  lNbVariables  =  3;
       //Array for the number of bins in each dimension :
     Int_t lNbBinsPerVar[3] = {0};
     lNbBinsPerVar[0] = 100;
     lNbBinsPerVar[1] = 1000;
     lNbBinsPerVar[2] = 22;
     fCFContCascadePIDAsOmegaMinus = new AliCFContainer(cfcontname_cascpidomegaminus,"Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
       //Setting the bin limits 
     fCFContCascadePIDAsOmegaMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
     fCFContCascadePIDAsOmegaMinus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
     fCFContCascadePIDAsOmegaMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
       //Setting the step title : one per PID case
     fCFContCascadePIDAsOmegaMinus->SetStepTitle(0, "No PID");
     fCFContCascadePIDAsOmegaMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
     fCFContCascadePIDAsOmegaMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
     fCFContCascadePIDAsOmegaMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
     fCFContCascadePIDAsOmegaMinus->SetStepTitle(4, "Comb. PID / Bachelor");
     fCFContCascadePIDAsOmegaMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
     fCFContCascadePIDAsOmegaMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
       //Setting the variable title, per axis
     fCFContCascadePIDAsOmegaMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
     fCFContCascadePIDAsOmegaMinus->SetVarTitle(1, "M( #Lambda , K^{-} ) (GeV/c^{2})");
     fCFContCascadePIDAsOmegaMinus->SetVarTitle(2, "Y_{#Omega}");
  }
  // -- PID container Omega+
  if(! fCFContCascadePIDAsOmegaPlus)  {
     const Int_t  lNbSteps      =  7;
     const Int_t  lNbVariables  =  3;
       //Array for the number of bins in each dimension :
     Int_t lNbBinsPerVar[3]= {0};
     lNbBinsPerVar[0] = 100;
     lNbBinsPerVar[1] = 1000;
     lNbBinsPerVar[2] = 22;  
     fCFContCascadePIDAsOmegaPlus = new AliCFContainer(cfcontname_cascpidomegaplus,"Pt_{cascade} Vs M_{#bar{#Omega}^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
       //Setting the bin limits 
     fCFContCascadePIDAsOmegaPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
     fCFContCascadePIDAsOmegaPlus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
     fCFContCascadePIDAsOmegaPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
       //Setting the step title : one per PID case
     fCFContCascadePIDAsOmegaPlus->SetStepTitle(0, "No PID");
     fCFContCascadePIDAsOmegaPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
     fCFContCascadePIDAsOmegaPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
     fCFContCascadePIDAsOmegaPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
     fCFContCascadePIDAsOmegaPlus->SetStepTitle(4, "Comb. PID / Bachelor");
     fCFContCascadePIDAsOmegaPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
     fCFContCascadePIDAsOmegaPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
       //Setting the variable title, per axis
     fCFContCascadePIDAsOmegaPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
     fCFContCascadePIDAsOmegaPlus->SetVarTitle(1, "M( #Lambda , K^{+} ) (GeV/c^{2})");
     fCFContCascadePIDAsOmegaPlus->SetVarTitle(2, "Y_{#Omega}");
  }
  // Container for optimisation of topological selections 
  if(! fCFContAsCascadeCuts){
	// Container meant to store all the relevant distributions corresponding to the cut variables.
        //          - NB overflow/underflow of variables on which we want to cut later should be 0!!!
     const Int_t  lNbSteps      =  4;
     const Int_t  lNbVariables  =  20;
       //Array for the number of bins in each dimension :
     Int_t lNbBinsPerVar[lNbVariables] = {0};
     lNbBinsPerVar[0]  = 25;   //DcaCascDaughters                : [0.0,2.,3.0]        -> Rec.Cut = 2.0; 
     lNbBinsPerVar[1]  = 25;   //DcaBachToPrimVertex             : [0.0,0.24,100.0]    -> Rec.Cur = 0.01;
     lNbBinsPerVar[2]  = 61;     //CascCosineOfPointingAngle    :  [0.94,1.01]          -> Rec.Cut = 0.95;
     //lNbBinsPerVar[2]  = 30;   //CascCosineOfPointingAngle       : [0.97,1.]           -> Rec.Cut = 0.98;
     lNbBinsPerVar[3]  = 40;   //CascRadius                      : [0.0,3.9,1000.0]    -> Rec.Cut = 0.2;
     lNbBinsPerVar[4]  = 30;   //InvMassLambdaAsCascDghter       : [1.1,1.3]           -> Rec.Cut = 0.008;
     lNbBinsPerVar[5]  = 20;   //DcaV0Daughters                  : [0.0,2.0]           -> Rec.Cut = 1.5;
     lNbBinsPerVar[6]  = 201;  //V0CosineOfPointingAngle         : [0.89,1.0]          -> Rec.Cut = 0.9;
     lNbBinsPerVar[7]  = 40;   //V0Radius                        : [0.0,3.9,1000.0]    -> Rec.Cut = 0.2;
     lNbBinsPerVar[8]  = 40;   //DcaV0ToPrimVertex               : [0.0,0.39,110.0]    -> Rec.Cut = 0.01;
     lNbBinsPerVar[9]  = 25;   //DcaPosToPrimVertex              : [0.0,0.24,100.0]    -> Rec.Cut = 0.05;
     lNbBinsPerVar[10] = 25;   //DcaNegToPrimVertex              : [0.0,0.24,100.0]    -> Rec.Cut = 0.05;
     lNbBinsPerVar[11] = 150;  //InvMassXi                       :  2-MeV/c2 bins
     lNbBinsPerVar[12] = 120;  //InvMassOmega                    :  2-MeV/c2 bins
     lNbBinsPerVar[13] = 100;  //CascTransvMom                   : [0.0,10.0]
     lNbBinsPerVar[14] = 110;  //Y(Xi)                           :  0.02 unit of y per bin 
     lNbBinsPerVar[15] = 110;  //Y(Omega)                        :  0.02 unit of y per bin
     lNbBinsPerVar[16] = 112;  //Proper lenght of cascade
     lNbBinsPerVar[17] = 112;  //Proper lenght of V0 
     lNbBinsPerVar[18] = 112;  //Distance V0-Xi in the transverse plane  
     lNbBinsPerVar[19] = 26;   //Proton cascade daughter momentum
     fCFContAsCascadeCuts = new AliCFContainer(cfcontname_casccuts,"Cut Container for Asso. Cascades", lNbSteps, lNbVariables, lNbBinsPerVar );
       //Setting the bin limits 
       //0 - DcaCascDaughters
     //Double_t *lBinLim0 = new Double_t[ lNbBinsPerVar[0]+1 ];
     //for(Int_t i=0; i<lNbBinsPerVar[0]; i++) lBinLim0[i] = (Double_t)0.0 + (2.4 -0.0)/(lNbBinsPerVar[0] - 1) * (Double_t)i;
     //lBinLim0[ lNbBinsPerVar[0] ] = 3.0;
     //fCFContAsCascadeCuts -> SetBinLimits(0, lBinLim0);
     //delete[] lBinLim0;
     fCFContAsCascadeCuts->SetBinLimits(0,0.0,2.5);
       //1 - DcaBachToPrimVertex
     Double_t *lBinLim1 = new Double_t[ lNbBinsPerVar[1]+1 ];
     for(Int_t i=0; i<lNbBinsPerVar[1]; i++) lBinLim1[i] = (Double_t)0.0 + (0.24 - 0.0)/(lNbBinsPerVar[1] - 1) * (Double_t)i;
     lBinLim1[ lNbBinsPerVar[1] ] = 100.0;
     fCFContAsCascadeCuts -> SetBinLimits(1, lBinLim1);
     delete [] lBinLim1;
       //2 - CascCosineOfPointingAngle
     fCFContAsCascadeCuts -> SetBinLimits(2, .94, 1.001);        
       //3 - CascRadius
     Double_t *lBinLim3 = new Double_t[ lNbBinsPerVar[3]+1 ];
     for(Int_t i=0; i<lNbBinsPerVar[3]; i++) lBinLim3[i] = (Double_t)0.0 + (3.9 -0.0)/(lNbBinsPerVar[3] - 1) * (Double_t)i;
     lBinLim3[ lNbBinsPerVar[3] ] = 1000.0;
     fCFContAsCascadeCuts -> SetBinLimits(3, lBinLim3);
     delete[] lBinLim3;
       //4 - InvMassLambdaAsCascDghter
     fCFContAsCascadeCuts->SetBinLimits(4, 1.1, 1.13); 
       //5 - DcaV0Daughters
     fCFContAsCascadeCuts->SetBinLimits(5, 0., 2.);        
       //6 - V0CosineOfPointingAngle
     fCFContAsCascadeCuts->SetBinLimits(6, 0.8, 1.001);
       //7 - V0Radius
     Double_t *lBinLim7 = new Double_t[ lNbBinsPerVar[7]+1 ];
     for(Int_t i=0; i<lNbBinsPerVar[7]; i++) lBinLim7[i] = (Double_t)0.0 + (3.9 - 0.0)/(lNbBinsPerVar[7] - 1) * (Double_t)i ;
     lBinLim7[ lNbBinsPerVar[7] ] = 1000.0;
     fCFContAsCascadeCuts -> SetBinLimits(7, lBinLim7);
     delete [] lBinLim7;      
       //8 - DcaV0ToPrimVertexXi : 0. to 0.4 
     Double_t *lBinLim8 = new Double_t[ lNbBinsPerVar[8]+1 ];
     for(Int_t i=0; i<lNbBinsPerVar[8]; i++) lBinLim8[i] = (Double_t)0.0 + (0.39 - 0.0)/(lNbBinsPerVar[8] - 1) * (Double_t)i ;
     lBinLim8[ lNbBinsPerVar[8] ] = 100.0;
     fCFContAsCascadeCuts -> SetBinLimits(8, lBinLim8);
     delete [] lBinLim8;      
       //9 - DcaPosToPrimVertexXi
     Double_t *lBinLim9 = new Double_t[ lNbBinsPerVar[9]+1 ];
     for(Int_t i=0; i<lNbBinsPerVar[9]; i++) lBinLim9[i] = (Double_t)0.0 + (0.24 - 0.0)/(lNbBinsPerVar[9] - 1) * (Double_t)i;
     lBinLim9[ lNbBinsPerVar[9] ] = 100.0;
     fCFContAsCascadeCuts -> SetBinLimits(9, lBinLim9);
     delete [] lBinLim9;   
       //10 - DcaNegToPrimVertexXi
     Double_t *lBinLim10 = new Double_t[ lNbBinsPerVar[10]+1 ];
     for(Int_t i=0; i<lNbBinsPerVar[10]; i++) lBinLim10[i] = (Double_t)0.0 + (0.24 - 0.0 )/(lNbBinsPerVar[10] - 1) * (Double_t)i;
     lBinLim10[ lNbBinsPerVar[10] ] = 100.0;
     fCFContAsCascadeCuts -> SetBinLimits(10, lBinLim10);
     delete [] lBinLim10;  
       //11 - InvMassXi
     fCFContAsCascadeCuts -> SetBinLimits(11, 1.25, 1.40);
       //12 - InvMassOmega
     fCFContAsCascadeCuts -> SetBinLimits(12, 1.62, 1.74);
       //13 - XiTransvMom 
     fCFContAsCascadeCuts -> SetBinLimits(13, 0.0, 10.0);
       //14 - Y(Xi) 
     fCFContAsCascadeCuts -> SetBinLimits(14, -1.1, 1.1);
       //15 - Y(Omega)
     fCFContAsCascadeCuts -> SetBinLimits(15, -1.1, 1.1); 
       //16 - Proper time cascade 
     Double_t *lBinLim16 = new Double_t[ lNbBinsPerVar[16]+1 ];
     for(Int_t i=0; i<lNbBinsPerVar[16]; i++) lBinLim16[i] = (Double_t)-1. + (110. + 1.0 )/(lNbBinsPerVar[16] - 1) * (Double_t)i;
     lBinLim16[ lNbBinsPerVar[16] ] = 2000.0;
     fCFContAsCascadeCuts -> SetBinLimits(16, lBinLim16);
       //17 - Proper time V0 
     fCFContAsCascadeCuts -> SetBinLimits(17, lBinLim16);
       //18 - Distance V0-Xi in the transverse plane
     fCFContAsCascadeCuts -> SetBinLimits(18, lBinLim16);
     delete [] lBinLim16;
      //19 - (Anti-)Proton cascade daughter momentum
     Double_t *lBinLim19  = new Double_t[ lNbBinsPerVar[19]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[19];i++)   lBinLim19[i] = (Double_t)0.0 + (0.5-0.0)/(lNbBinsPerVar[19]-1) * (Double_t)i;
        lBinLim19[ lNbBinsPerVar[19]  ] = 50.0;
     fCFContAsCascadeCuts->SetBinLimits(19, lBinLim19);
     delete [] lBinLim19;
       // Setting the number of steps : one for each cascade species (Xi-, Xi+ and Omega-, Omega+)
     fCFContAsCascadeCuts->SetStepTitle(0, "#Xi^{-} candidates associated to MC");
     fCFContAsCascadeCuts->SetStepTitle(1, "#bar{#Xi}^{+} candidates associated to MC");
     fCFContAsCascadeCuts->SetStepTitle(2, "#Omega^{-} candidates associated to MC");
     fCFContAsCascadeCuts->SetStepTitle(3, "#bar{#Omega}^{+} candidates associated to MC");
       // Setting the variable title, per axis
     fCFContAsCascadeCuts->SetVarTitle(0,  "DCA(cascade daughters) (cm)");
     fCFContAsCascadeCuts->SetVarTitle(1,  "ImpactParamToPV(bachelor) (cm)");
     fCFContAsCascadeCuts->SetVarTitle(2,  "cos(cascade PA)");
     fCFContAsCascadeCuts->SetVarTitle(3,  "R_{2d}(cascade decay) (cm)");
     fCFContAsCascadeCuts->SetVarTitle(4,  "M_{#Lambda}(as casc dghter) (GeV/c^{2})");
     fCFContAsCascadeCuts->SetVarTitle(5,  "DCA(V0 daughters) (cm)");
     if      (fCollidingSystem == 0) fCFContAsCascadeCuts->SetVarTitle(6,  "cos(V0 PA) to cascade vtx");
     else if (fCollidingSystem == 1) fCFContAsCascadeCuts->SetVarTitle(6,  "cos(V0 PA) to primary vtx");
     fCFContAsCascadeCuts->SetVarTitle(7,  "R_{2d}(V0 decay) (cm)");
     fCFContAsCascadeCuts->SetVarTitle(8,  "ImpactParamToPV(V0) (cm)");
     fCFContAsCascadeCuts->SetVarTitle(9,  "ImpactParamToPV(Pos) (cm)");
     fCFContAsCascadeCuts->SetVarTitle(10, "ImpactParamToPV(Neg) (cm)");
     fCFContAsCascadeCuts->SetVarTitle(11, "Inv. Mass(Xi) (GeV/c^{2})");
     fCFContAsCascadeCuts->SetVarTitle(12, "Inv. Mass(Omega) (GeV/c^{2})");
     fCFContAsCascadeCuts->SetVarTitle(13, "Pt_{MC}(cascade) (GeV/c)");
     fCFContAsCascadeCuts->SetVarTitle(14, "Y_{MC}(Xi)");
     fCFContAsCascadeCuts->SetVarTitle(15, "Y_{MC}(Omega)");
     fCFContAsCascadeCuts->SetVarTitle(16, "mL/p cascade (cm)");
     fCFContAsCascadeCuts->SetVarTitle(17, "mL/p V0 (cm)"); 
     fCFContAsCascadeCuts->SetVarTitle(18, "Distance V0-Cascade in the transverse plane (cm)");
     fCFContAsCascadeCuts->SetVarTitle(19, "(Anti-)Proton casc. daught. momemtum (GeV/c)");
  }

 PostData(1, fListHistCascade); 
 PostData(2, fCFContCascadePIDAsXiMinus);
 PostData(3, fCFContCascadePIDAsXiPlus);
 PostData(4, fCFContCascadePIDAsOmegaMinus);
 PostData(5, fCFContCascadePIDAsOmegaPlus);
 PostData(6, fCFContAsCascadeCuts);

}// end CreateOutputObjects


//________________________________________________________________________
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// - EXECUTE ALL THE CALCULATION AND FILL THE OUTPUTS
void AliAnalysisTaskCheckPerformanceCascadepp::UserExec(Option_t *) {
	
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Main loop (called for each event)
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   //------------------
   // - Define variables
   AliESDEvent  *lESDevent = 0x0;
   AliAODEvent  *lAODevent = 0x0;
   AliMCEvent   *lMCevent  = 0x0; 
   AliStack     *lMCstack  = 0x0; 
   TClonesArray *arrayMC = 0;

   //-------------------------
   // - Check the PID response
   if (!fPIDResponse) {
        AliError("Cannot get pid response");
        return;
   }

  ///////////////////
  // EVENT SELECTIONS
  ///////////////////
  // In order:
  // 0) SDD status 
  //       <---------- (A) Check generated particle pT spectra before all the event selections
  // 1) Pre-Trigger selections
  //     1.1) Incomplete DAQ events (introduced for Run2 2015 data)
  //     1.2) Background rejection based on SPD cluster vs tracklet correlation
  //     1.3) Pileup
  //       <---------- (B) Fill efficiency denominator + check generated particle pT spectra before the physics selection
  // 2) Trigger selection (Physics selection)
  //       <---------- (C) Check generated particle pT spectra  after the physics selection
  // -) Cascade and V0 re-vertexer
  //       <---------- (D) Check generated particle pT spectra after the re-vertexing
  // 3) Well-established PV
  //     3.1) not only TPC vertex
  //       <---------- (E) Check generated particle pT spectra after the request on both SPD and TPC vertices
  //     3.2) requirement on the resolution and dispersion
  //       <---------- (F) Check generated particle pT spectra after the requests on resolution and dispersion
  //     3.2) distance between the two vertices
  //       <---------- (G) Check generated particle pT spectra after the vertices proximity check
  //     3.4) |Zpv| < 10 cm
  //    <---------- (H) Fill the efficiency denominator after all the event selections
  // - Define useful variables
  Int_t ncascades          = 0;
  Int_t nTrackMultiplicity = 0;
    
   //=======================================================
   // Load the InputEvent and check it (for the ESD and AOD)
   //=======================================================
   if (fAnalysisType == "ESD") {
       lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
       if (!lESDevent) {
	   Printf("ERROR: lESDevent not available \n");
	   cout << "Name of the file with pb :" <<  CurrentFileName() << endl;  
	   return;
       }
       lMCevent = MCEvent();
       if (!lMCevent) {
           Printf("ERROR: Could not retrieve MC event \n");
           cout << "Name of the file with pb :" <<  CurrentFileName() << endl;
           return;
       }
       lMCstack = lMCevent->Stack();
       if (!lMCstack) {
           Printf("ERROR: Could not retrieve MC stack \n");
           cout << "Name of the file with pb :" <<  CurrentFileName() << endl;
           return;
       }
   } else if (fAnalysisType == "AOD") {  
       lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() ); 
       if (!lAODevent) {
       	   Printf("ERROR: lAODevent not available \n");
	   cout << "Name of the file with pb :" <<  CurrentFileName() << endl;
	   return;
       }
       arrayMC = (TClonesArray*) lAODevent->GetList()->FindObject(AliAODMCParticle::StdBranchName());
       if (!arrayMC) AliFatal("Error: MC particles branch not found!\n");
   } else {
       Printf("Analysis type (ESD or AOD) not specified \n");
       return;
   }

   //============================
   // Plots Before any selections
   //============================
   if (fAnalysisType == "ESD") {
       nTrackMultiplicity = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5);
       ncascades = lESDevent->GetNumberOfCascades();
   } else if (fAnalysisType == "AOD") {
       nTrackMultiplicity = -100;
       ncascades = lAODevent->GetNumberOfCascades();
   }
   fHistCascadeMultiplicityBeforeAnySel->Fill(ncascades);
   

   //========================
   // 0) SDD status selection
   //========================
   if (fApplyEvSelSDDstatus && fCollidingSystem == 0) {
        TString trcl = " ";
        if      (fAnalysisType == "ESD") trcl = lESDevent->GetFiredTriggerClasses();
        else if (fAnalysisType == "AOD") trcl = lAODevent->GetFiredTriggerClasses();
        if (fwithSDD && !(trcl.Contains("ALLNOTRD"))) {
                 AliWarning("We are selecting events with SDD turn ON. This event has the SDD turn OFF. =>  RETURN!! (Exclude it)...");
                 PostData(1, fListHistCascade);
                 PostData(2, fCFContCascadePIDAsXiMinus);
                 PostData(3, fCFContCascadePIDAsXiPlus);
                 PostData(4, fCFContCascadePIDAsOmegaMinus);
                 PostData(5, fCFContCascadePIDAsOmegaPlus);
                 PostData(6, fCFContAsCascadeCuts);
                 return;
        } else if (!fwithSDD && (trcl.Contains("ALLNOTRD"))) {
                 AliWarning("We are selecting events with SDD turn OFF. This event has the SDD turn ON. =>  RETURN!! (Exclude it)...");
                 PostData(1, fListHistCascade);
                 PostData(2, fCFContCascadePIDAsXiMinus);
                 PostData(3, fCFContCascadePIDAsXiPlus);
                 PostData(4, fCFContCascadePIDAsOmegaMinus);
                 PostData(5, fCFContCascadePIDAsOmegaPlus);
                 PostData(6, fCFContAsCascadeCuts);
                 return;
        }
   }
   fHistCascadeMultiplicityAfterSDDstatusSel->Fill(ncascades);

   // - (A) Check generated particle pT spectra before all the event selections
   Int_t lNbMCPrimary_A = 0;   lNbMCPrimary_A = lMCstack->GetNprimary();
   Int_t lPrimaryTrackMultiplicity_A = nTrackMultiplicity;
   Double_t partEnergy_A, partPz_A, partPt_A, lRapXiMC_A;
   Int_t PDGcode_A;
   TParticle* lCurrentParticlePrimary_A = 0x0;        
   AliAODMCParticle *lCurrentParticleaod_A = 0x0;
   for (Int_t iCurrentLabelStack_A = 0; iCurrentLabelStack_A < lNbMCPrimary_A; iCurrentLabelStack_A++) {
         partEnergy_A = 0.;  partPz_A = 0.;  partPt_A  = 0.;  PDGcode_A = 0;
         if (fAnalysisType == "ESD") {
              lCurrentParticlePrimary_A = 0x0;  lCurrentParticlePrimary_A = lMCstack->Particle(iCurrentLabelStack_A);        
              if (!lCurrentParticlePrimary_A) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack_A)) continue;
              partEnergy_A = lCurrentParticlePrimary_A->Energy();   partPz_A = lCurrentParticlePrimary_A->Pz();  partPt_A = lCurrentParticlePrimary_A->Pt();  PDGcode_A = lCurrentParticlePrimary_A->GetPdgCode();
         } else if (fAnalysisType == "AOD") {
              lCurrentParticleaod_A = 0x0;  lCurrentParticleaod_A = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack_A);
              if (!lCurrentParticleaod_A) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lCurrentParticleaod_A->IsPhysicalPrimary()) continue;
              partEnergy_A = lCurrentParticleaod_A->E();  partPz_A = lCurrentParticleaod_A->Pz();  partPt_A = lCurrentParticleaod_A->Pt();  PDGcode_A = lCurrentParticleaod_A->GetPdgCode();
         }
         lRapXiMC_A = 0.5*TMath::Log((partEnergy_A + partPz_A) / (partEnergy_A - partPz_A + 1.e-13));
         if (PDGcode_A ==  3312)  f3dHistGenPtVsGenYvsNtracksXiMinus_A->Fill(partPt_A, lRapXiMC_A, lPrimaryTrackMultiplicity_A);
         if (PDGcode_A == -3312)  f3dHistGenPtVsGenYvsNtracksXiPlus_A->Fill(partPt_A, lRapXiMC_A, lPrimaryTrackMultiplicity_A);
         if (PDGcode_A ==  3334)  f3dHistGenPtVsGenYvsNtracksOmegaMinus_A->Fill(partPt_A, lRapXiMC_A, lPrimaryTrackMultiplicity_A);
         if (PDGcode_A == -3334)  f3dHistGenPtVsGenYvsNtracksOmegaPlus_A->Fill(partPt_A, lRapXiMC_A, lPrimaryTrackMultiplicity_A);
   }

   //======================================
   // 1.1) Removal of incomplete DAQ events
   //======================================
   // - Incomplete DAQ events rejection Run2 data 2015
   Bool_t IncompleteDAQ = kFALSE;
   if (fAnalysisType == "ESD") IncompleteDAQ = lESDevent->IsIncompleteDAQ();
   if (fApplyEvSelDAQincomplete && IncompleteDAQ){
        AliWarning("This is a DAQ incomplete event... return !");
        PostData(1, fListHistCascade);
        PostData(2, fCFContCascadePIDAsXiMinus);
        PostData(3, fCFContCascadePIDAsXiPlus);
        PostData(4, fCFContCascadePIDAsOmegaMinus);
        PostData(5, fCFContCascadePIDAsOmegaPlus);
        PostData(6, fCFContAsCascadeCuts);
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
              PostData(2, fCFContCascadePIDAsXiMinus);
              PostData(3, fCFContCascadePIDAsXiPlus);
              PostData(4, fCFContCascadePIDAsOmegaMinus);
              PostData(5, fCFContCascadePIDAsOmegaPlus);
              PostData(6, fCFContAsCascadeCuts);
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
              PostData(2, fCFContCascadePIDAsXiMinus);
              PostData(3, fCFContCascadePIDAsXiPlus);
              PostData(4, fCFContCascadePIDAsOmegaMinus);
              PostData(5, fCFContCascadePIDAsOmegaPlus);
              PostData(6, fCFContAsCascadeCuts);
              return;
           }
      } else if (fAnalysisType == "AOD") {
           if(lAODevent->IsPileupFromSPD(fSPDPileUpminContributors)){
              AliWarning("Pb / Pile-up event ... return!");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDAsXiMinus);
              PostData(3, fCFContCascadePIDAsXiPlus);
              PostData(4, fCFContCascadePIDAsOmegaMinus);
              PostData(5, fCFContCascadePIDAsOmegaPlus);
              PostData(6, fCFContAsCascadeCuts);
              return;
           }
      }
   }
   fHistCascadeMultiplicityAfterPileupRej->Fill(ncascades);

   // - (B) Efficiency denominator before physics selection 
   Int_t lNbMCPrimary_B = 0;   lNbMCPrimary_B = lMCstack->GetNprimary();
   Int_t lPrimaryTrackMultiplicity_B = nTrackMultiplicity;
   Double_t partEnergy_B, partPz_B, partP_B, partPt_B, partVx_B, partVy_B, partVz_B, bacVx_B, bacVy_B, bacVz_B, partMass_B, lRapXiMC_B, lctau_B;
   Int_t PDGcode_B;
   TParticle* lCurrentParticlePrimary_B = 0x0;     TParticle *mcBach_B = 0x0;
   AliAODMCParticle *lCurrentParticleaod_B = 0x0;  AliAODMCParticle *mcBachaod_B = 0x0;
   for (Int_t iCurrentLabelStack_B = 0; iCurrentLabelStack_B < lNbMCPrimary_B; iCurrentLabelStack_B++) {
         partEnergy_B = 0.;  partPz_B = 0.;  partP_B    = 0.;  partPt_B  = 0.;
         partVx_B     = 0.;  partVy_B = 0.;  partVz_B   = 0.;  bacVx_B   = 0.;
         bacVy_B      = 0.;  bacVz_B  = 0.;  partMass_B = 0.;  PDGcode_B = 0;
         if (fAnalysisType == "ESD") {
              lCurrentParticlePrimary_B = 0x0;  lCurrentParticlePrimary_B = lMCstack->Particle(iCurrentLabelStack_B);
              if (!lCurrentParticlePrimary_B) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack_B)) continue;
              partEnergy_B = lCurrentParticlePrimary_B->Energy();   partPz_B = lCurrentParticlePrimary_B->Pz();  partPt_B = lCurrentParticlePrimary_B->Pt();  partP_B  = lCurrentParticlePrimary_B->P();
              partMass_B   = lCurrentParticlePrimary_B->GetMass();  partVx_B = lCurrentParticlePrimary_B->Vx();  partVy_B = lCurrentParticlePrimary_B->Vy();  partVz_B = lCurrentParticlePrimary_B->Vz();
              PDGcode_B = lCurrentParticlePrimary_B->GetPdgCode();
              if (lCurrentParticlePrimary_B->GetDaughter(0) >= 0) {
                  mcBach_B = 0x0;  mcBach_B = lMCstack->Particle(lCurrentParticlePrimary_B->GetDaughter(0));
                  if (mcBach_B) { bacVx_B = mcBach_B->Vx();  bacVy_B = mcBach_B->Vy();  bacVz_B = mcBach_B->Vz(); }
              }
         } else if (fAnalysisType == "AOD") {
              lCurrentParticleaod_B = 0x0;  lCurrentParticleaod_B = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack_B);
              if (!lCurrentParticleaod_B) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lCurrentParticleaod_B->IsPhysicalPrimary()) continue;
              partEnergy_B = lCurrentParticleaod_B->E();  partPz_B = lCurrentParticleaod_B->Pz();  partP_B  = lCurrentParticleaod_B->P();   partPt_B = lCurrentParticleaod_B->Pt();
              partMass_B   = lCurrentParticleaod_B->M();  partVx_B = lCurrentParticleaod_B->Xv();  partVy_B = lCurrentParticleaod_B->Yv();  partVz_B = lCurrentParticleaod_B->Zv();
              PDGcode_B = lCurrentParticleaod_B->GetPdgCode();
              if (lCurrentParticleaod_B->GetDaughter(0) >= 0) {
                   mcBachaod_B = 0x0;  mcBachaod_B = (AliAODMCParticle*) arrayMC->At(lCurrentParticleaod_B->GetDaughter(0));
                   if (mcBachaod_B) { bacVx_B = mcBachaod_B->Xv();  bacVy_B = mcBachaod_B->Yv();  bacVz_B = mcBachaod_B->Zv(); }
              } 
         }
         lRapXiMC_B = 0.5*TMath::Log((partEnergy_B + partPz_B) / (partEnergy_B - partPz_B + 1.e-13));
         lctau_B    = TMath::Sqrt((partVx_B-bacVx_B)*(partVx_B-bacVx_B)+(partVy_B-bacVy_B)*(partVy_B-bacVy_B)+(partVz_B-bacVz_B)*(partVz_B-bacVz_B));
         if (partP_B != 0.) lctau_B = lctau_B*partMass_B/partP_B;
         else               lctau_B = -1.;
         if (PDGcode_B ==  3312) {
             f3dHistGenPtVsGenYvsNtracksXiMinus_B->Fill(partPt_B, lRapXiMC_B, lPrimaryTrackMultiplicity_B);
             f3dHistGenPtVsGenctauvsYXiMinus_B->Fill(partPt_B, lctau_B, lRapXiMC_B);
         }
         if (PDGcode_B == -3312) {
             f3dHistGenPtVsGenYvsNtracksXiPlus_B->Fill(partPt_B, lRapXiMC_B, lPrimaryTrackMultiplicity_B);
             f3dHistGenPtVsGenctauvsYXiPlus_B->Fill(partPt_B, lctau_B, lRapXiMC_B);
         }
         if (PDGcode_B ==  3334) {
             f3dHistGenPtVsGenYvsNtracksOmegaMinus_B->Fill(partPt_B, lRapXiMC_B, lPrimaryTrackMultiplicity_B);
             f3dHistGenPtVsGenctauvsYOmegaMinus_B->Fill(partPt_B, lctau_B, lRapXiMC_B);
         }
         if (PDGcode_B == -3334) {
             f3dHistGenPtVsGenYvsNtracksOmegaPlus_B->Fill(partPt_B, lRapXiMC_B, lPrimaryTrackMultiplicity_B);
             f3dHistGenPtVsGenctauvsYOmegaPlus_B->Fill(partPt_B, lctau_B, lRapXiMC_B);
         }
   }

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
           PostData(2, fCFContCascadePIDAsXiMinus);
           PostData(3, fCFContCascadePIDAsXiPlus);
           PostData(4, fCFContCascadePIDAsOmegaMinus);
           PostData(5, fCFContCascadePIDAsOmegaPlus);
           PostData(6, fCFContAsCascadeCuts);
           return;
       }
   }
   fHistCascadeMultiplicityAfterPhysicsSel->Fill(ncascades);

   // - (C) Efficiency denominator after physics selection 
   Int_t lNbMCPrimary_C = 0;   lNbMCPrimary_C = lMCstack->GetNprimary();
   Int_t lPrimaryTrackMultiplicity_C = nTrackMultiplicity;
   Double_t partEnergy_C, partPz_C, partPt_C, lRapXiMC_C;
   Int_t PDGcode_C;
   TParticle* lCurrentParticlePrimary_C = 0x0;                               
   AliAODMCParticle *lCurrentParticleaod_C = 0x0;  
   for (Int_t iCurrentLabelStack_C = 0; iCurrentLabelStack_C < lNbMCPrimary_C; iCurrentLabelStack_C++) {
         partEnergy_C = 0.;  partPz_C = 0.;  partPt_C = 0.;  PDGcode_C = 0;
         if (fAnalysisType == "ESD") {
              lCurrentParticlePrimary_C = 0x0;  lCurrentParticlePrimary_C = lMCstack->Particle(iCurrentLabelStack_C);
              if (!lCurrentParticlePrimary_C) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack_C)) continue;
              partEnergy_C = lCurrentParticlePrimary_C->Energy();  partPz_C = lCurrentParticlePrimary_C->Pz();  partPt_C = lCurrentParticlePrimary_C->Pt();  PDGcode_C = lCurrentParticlePrimary_C->GetPdgCode();
         } else if (fAnalysisType == "AOD") {
              lCurrentParticleaod_C = 0x0;  lCurrentParticleaod_C = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack_C);
              if (!lCurrentParticleaod_C) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lCurrentParticleaod_C->IsPhysicalPrimary()) continue;
              partEnergy_C = lCurrentParticleaod_C->E();  partPz_C = lCurrentParticleaod_C->Pz();  partPt_C = lCurrentParticleaod_C->Pt();  PDGcode_C = lCurrentParticleaod_C->GetPdgCode();
         }
         lRapXiMC_C = 0.5*TMath::Log((partEnergy_C + partPz_C) / (partEnergy_C - partPz_C + 1.e-13));
         if (PDGcode_C ==  3312)  f3dHistGenPtVsGenYvsNtracksXiMinus_C->Fill(partPt_C, lRapXiMC_C, lPrimaryTrackMultiplicity_C);
         if (PDGcode_C == -3312)  f3dHistGenPtVsGenYvsNtracksXiPlus_C->Fill(partPt_C, lRapXiMC_C, lPrimaryTrackMultiplicity_C);
         if (PDGcode_C ==  3334)  f3dHistGenPtVsGenYvsNtracksOmegaMinus_C->Fill(partPt_C, lRapXiMC_C, lPrimaryTrackMultiplicity_C);
         if (PDGcode_C == -3334)  f3dHistGenPtVsGenYvsNtracksOmegaPlus_C->Fill(partPt_C, lRapXiMC_C, lPrimaryTrackMultiplicity_C);
   }

   //===============================================
   // - Re-run V0 and cascade vertexers (only for ESD)
   //===============================================
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

   // - (D) Efficiency denominator after Re-vertexing
   Int_t lNbMCPrimary_D = 0;   lNbMCPrimary_D = lMCstack->GetNprimary();
   Int_t lPrimaryTrackMultiplicity_D = nTrackMultiplicity;
   Double_t partEnergy_D, partPz_D, partPt_D, lRapXiMC_D;
   Int_t PDGcode_D;
   TParticle* lCurrentParticlePrimary_D = 0x0;                               
   AliAODMCParticle *lCurrentParticleaod_D = 0x0;  
   for (Int_t iCurrentLabelStack_D = 0; iCurrentLabelStack_D < lNbMCPrimary_D; iCurrentLabelStack_D++) {
         partEnergy_D = 0.;  partPz_D = 0.;  partPt_D  = 0.;  PDGcode_D = 0;
         if (fAnalysisType == "ESD") {
              lCurrentParticlePrimary_D = 0x0;  lCurrentParticlePrimary_D = lMCstack->Particle(iCurrentLabelStack_D);
              if (!lCurrentParticlePrimary_D) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack_D)) continue;
              partEnergy_D = lCurrentParticlePrimary_D->Energy();  partPz_D = lCurrentParticlePrimary_D->Pz();  partPt_D = lCurrentParticlePrimary_D->Pt();  PDGcode_D = lCurrentParticlePrimary_D->GetPdgCode();
         } else if (fAnalysisType == "AOD") {
              lCurrentParticleaod_D = 0x0;  lCurrentParticleaod_D = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack_D);
              if (!lCurrentParticleaod_D) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lCurrentParticleaod_D->IsPhysicalPrimary()) continue;
              partEnergy_D = lCurrentParticleaod_D->E();  partPz_D = lCurrentParticleaod_D->Pz();  partPt_D = lCurrentParticleaod_D->Pt();  PDGcode_D = lCurrentParticleaod_D->GetPdgCode();
         }
         lRapXiMC_D = 0.5*TMath::Log((partEnergy_D + partPz_D) / (partEnergy_D - partPz_D + 1.e-13));
         if (PDGcode_D ==  3312)  f3dHistGenPtVsGenYvsNtracksXiMinus_D->Fill(partPt_D, lRapXiMC_D, lPrimaryTrackMultiplicity_D);
         if (PDGcode_D == -3312)  f3dHistGenPtVsGenYvsNtracksXiPlus_D->Fill(partPt_D, lRapXiMC_D, lPrimaryTrackMultiplicity_D);
         if (PDGcode_D ==  3334)  f3dHistGenPtVsGenYvsNtracksOmegaMinus_D->Fill(partPt_D, lRapXiMC_D, lPrimaryTrackMultiplicity_D);
         if (PDGcode_D == -3334)  f3dHistGenPtVsGenYvsNtracksOmegaPlus_D->Fill(partPt_D, lRapXiMC_D, lPrimaryTrackMultiplicity_D);
   }

   //=================================
   // 3) Well-established PV selection
   //=================================
   // - Just take the vertices
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
 
   //===============================================================================================
   // 3.1) reject events if both SPD and TPC vertices are explicitly requested and none is available
   //===============================================================================================
   if (fCollidingSystem == 0) {    
       if (fApplyEvSelNoTPConlyPrimVtx) {
            if (!(lESDPrimarySPDVtx->GetStatus() && lESDPrimaryTrackingVtx->GetStatus()) && fAnalysisType == "ESD"){
                  AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDAsXiMinus);
                  PostData(3, fCFContCascadePIDAsXiPlus);
                  PostData(4, fCFContCascadePIDAsOmegaMinus);
                  PostData(5, fCFContCascadePIDAsOmegaPlus);
                  PostData(6, fCFContAsCascadeCuts);
                  return;
            }
            if (!(lAODPrimarySPDVtx && lAODPrimaryTrackingVtx) && fAnalysisType == "AOD") {
                  AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDAsXiMinus);
                  PostData(3, fCFContCascadePIDAsXiPlus);
                  PostData(4, fCFContCascadePIDAsOmegaMinus);
                  PostData(5, fCFContCascadePIDAsOmegaPlus);
                  PostData(6, fCFContAsCascadeCuts);
                  return;
             }
        }
   }
   fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel->Fill(ncascades);

   // - (E) Efficiency denominator after rejection of events that has only vertex from TPC
   Int_t lNbMCPrimary_E = 0;  lNbMCPrimary_E = lMCstack->GetNprimary();
   Int_t lPrimaryTrackMultiplicity_E = nTrackMultiplicity;
   Double_t partEnergy_E, partPz_E, partPt_E, lRapXiMC_E;
   Int_t PDGcode_E;
   TParticle* lCurrentParticlePrimary_E = 0x0;
   AliAODMCParticle *lCurrentParticleaod_E = 0x0;
   for (Int_t iCurrentLabelStack_E = 0; iCurrentLabelStack_E < lNbMCPrimary_E; iCurrentLabelStack_E++) {
         partEnergy_E = 0.;  partPz_E = 0.;  partPt_E = 0.;  PDGcode_E = 0;
         if (fAnalysisType == "ESD") {
              lCurrentParticlePrimary_E = 0x0;  lCurrentParticlePrimary_E = lMCstack->Particle(iCurrentLabelStack_E);
              if (!lCurrentParticlePrimary_E) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack_E)) continue;
              partEnergy_E = lCurrentParticlePrimary_E->Energy();  partPz_E = lCurrentParticlePrimary_E->Pz();  partPt_E = lCurrentParticlePrimary_E->Pt();  PDGcode_E = lCurrentParticlePrimary_E->GetPdgCode();
         } else if (fAnalysisType == "AOD") {
              lCurrentParticleaod_E = 0x0;  lCurrentParticleaod_E = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack_E);
              if (!lCurrentParticleaod_E) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lCurrentParticleaod_E->IsPhysicalPrimary()) continue;
              partEnergy_E = lCurrentParticleaod_E->E();  partPz_E = lCurrentParticleaod_E->Pz();  partPt_E = lCurrentParticleaod_E->Pt();  PDGcode_E = lCurrentParticleaod_E->GetPdgCode();
         }
         lRapXiMC_E = 0.5*TMath::Log((partEnergy_E + partPz_E) / (partEnergy_E - partPz_E + 1.e-13));
         if (PDGcode_E ==  3312)  f3dHistGenPtVsGenYvsNtracksXiMinus_E->Fill(partPt_E, lRapXiMC_E, lPrimaryTrackMultiplicity_E);
         if (PDGcode_E == -3312)  f3dHistGenPtVsGenYvsNtracksXiPlus_E->Fill(partPt_E, lRapXiMC_E, lPrimaryTrackMultiplicity_E);
         if (PDGcode_E ==  3334)  f3dHistGenPtVsGenYvsNtracksOmegaMinus_E->Fill(partPt_E, lRapXiMC_E, lPrimaryTrackMultiplicity_E);
         if (PDGcode_E == -3334)  f3dHistGenPtVsGenYvsNtracksOmegaPlus_E->Fill(partPt_E, lRapXiMC_E, lPrimaryTrackMultiplicity_E);
   }

   //=================================================================
   // 3.2) check the spd vertex resolution and reject if not satisfied //FIXME: only for ESD
   //=================================================================
   if (!lESDPrimaryTrackingVtx->GetStatus()) {
        if (!lESDPrimarySPDVtx->GetStatus()) {
             AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
             PostData(1, fListHistCascade);
             PostData(2, fCFContCascadePIDAsXiMinus);
             PostData(3, fCFContCascadePIDAsXiPlus);
             PostData(4, fCFContCascadePIDAsOmegaMinus);
             PostData(5, fCFContCascadePIDAsOmegaPlus);
             PostData(6, fCFContAsCascadeCuts);
             return;
        }
        if (fApplyEvSelSPDvtxres && lESDPrimarySPDVtx->GetStatus() && lESDPrimarySPDVtx->IsFromVertexerZ() && !(lESDPrimarySPDVtx->GetDispersion()<0.04 && lESDPrimarySPDVtx->GetZRes()<0.25)) {
             AliWarning("Pb / The SPD prim. vertex has a Z resolution > 0.25 and dispersion > 0.04 ... return !");
             PostData(1, fListHistCascade);
             PostData(2, fCFContCascadePIDAsXiMinus);
             PostData(3, fCFContCascadePIDAsXiPlus);
             PostData(4, fCFContCascadePIDAsOmegaMinus);
             PostData(5, fCFContCascadePIDAsOmegaPlus);
             PostData(6, fCFContAsCascadeCuts);
             return;
        }
   } else {
        if (fApplyEvSelSPDvtxres && lESDPrimarySPDVtx->GetStatus() && lESDPrimarySPDVtx->IsFromVertexerZ() && !(lESDPrimarySPDVtx->GetDispersion()<0.04 && lESDPrimarySPDVtx->GetZRes()<0.25)) {
             AliWarning("Pb / The SPD prim. vertex has a Z resolution > 0.25 and dispersion > 0.04 ... return !");
             PostData(1, fListHistCascade);
             PostData(2, fCFContCascadePIDAsXiMinus);
             PostData(3, fCFContCascadePIDAsXiPlus);
             PostData(4, fCFContCascadePIDAsOmegaMinus);
             PostData(5, fCFContCascadePIDAsOmegaPlus);
             PostData(6, fCFContAsCascadeCuts);
             return;
        }
   }
   fHistCascadeMultiplicityAfterSPDresolution->Fill(ncascades);

   // - (F) Efficiency denominator after the rejection of events that do not satisfy the request on the dispersion and resolution of the SPD vertex
   Int_t lNbMCPrimary_F = 0;   lNbMCPrimary_F = lMCstack->GetNprimary();
   Int_t lPrimaryTrackMultiplicity_F = nTrackMultiplicity;
   Double_t partEnergy_F, partPz_F, partPt_F, lRapXiMC_F;
   Int_t PDGcode_F;
   TParticle* lCurrentParticlePrimary_F = 0x0;
   AliAODMCParticle *lCurrentParticleaod_F = 0x0;
   for (Int_t iCurrentLabelStack_F = 0; iCurrentLabelStack_F < lNbMCPrimary_F; iCurrentLabelStack_F++) {
         partEnergy_F = 0.;  partPz_F = 0.;  partPt_F = 0.;  PDGcode_F = 0;
         if (fAnalysisType == "ESD") {
              lCurrentParticlePrimary_F = 0x0;  lCurrentParticlePrimary_F = lMCstack->Particle(iCurrentLabelStack_F);
              if (!lCurrentParticlePrimary_F) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack_F)) continue;
              partEnergy_F = lCurrentParticlePrimary_F->Energy();  partPz_F = lCurrentParticlePrimary_F->Pz();  partPt_F = lCurrentParticlePrimary_F->Pt();  PDGcode_F = lCurrentParticlePrimary_F->GetPdgCode();
         } else if (fAnalysisType == "AOD") {
              lCurrentParticleaod_F = 0x0;  lCurrentParticleaod_F = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack_F);
              if (!lCurrentParticleaod_F) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lCurrentParticleaod_F->IsPhysicalPrimary()) continue;
              partEnergy_F = lCurrentParticleaod_F->E();  partPz_F = lCurrentParticleaod_F->Pz();  partPt_F = lCurrentParticleaod_F->Pt();  PDGcode_F = lCurrentParticleaod_F->GetPdgCode();
         }
         lRapXiMC_F = 0.5*TMath::Log((partEnergy_F + partPz_F) / (partEnergy_F - partPz_F + 1.e-13));
         if (PDGcode_F ==  3312)  f3dHistGenPtVsGenYvsNtracksXiMinus_F->Fill(partPt_F, lRapXiMC_F, lPrimaryTrackMultiplicity_F);
         if (PDGcode_F == -3312)  f3dHistGenPtVsGenYvsNtracksXiPlus_F->Fill(partPt_F, lRapXiMC_F, lPrimaryTrackMultiplicity_F);
         if (PDGcode_F ==  3334)  f3dHistGenPtVsGenYvsNtracksOmegaMinus_F->Fill(partPt_F, lRapXiMC_F, lPrimaryTrackMultiplicity_F);
         if (PDGcode_F == -3334)  f3dHistGenPtVsGenYvsNtracksOmegaPlus_F->Fill(partPt_F, lRapXiMC_F, lPrimaryTrackMultiplicity_F);
   }

   //==============================================================================================
   // 3.3) check the proximity between the spd vertex and track vertex, and reject if not satisfied
   //============================================================================================== 
   if (fCollidingSystem == 0) {
        if (lESDPrimaryTrackingVtx->GetStatus() && lESDPrimarySPDVtx->GetStatus()) {
             if (fApplyEvSelVtxProximity && (TMath::Abs(lESDPrimarySPDVtx->GetZ() - lESDPrimaryTrackingVtx->GetZ()) > 0.5)) {
                  PostData(1, fListHistCascade);
                  PostData(2, fCFContCascadePIDAsXiMinus);
                  PostData(3, fCFContCascadePIDAsXiPlus);
                  PostData(4, fCFContCascadePIDAsOmegaMinus);
                  PostData(5, fCFContCascadePIDAsOmegaPlus);
                  PostData(6, fCFContAsCascadeCuts);
                  return;
             }
        }
   } 
   // - vertex selection for pPb analysis
   else if (fCollidingSystem == 1) {
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
                   PostData(2, fCFContCascadePIDAsXiMinus);
                   PostData(3, fCFContCascadePIDAsXiPlus);
                   PostData(4, fCFContCascadePIDAsOmegaMinus);
                   PostData(5, fCFContCascadePIDAsOmegaPlus);
                   PostData(6, fCFContAsCascadeCuts);
                   return;
              }
              if (fUtils->IsFirstEventInChunk(lESDevent)) { //Is First event in chunk rejection: Still present!
                   AliWarning("Pb / This is the first event in the chunk! ... return !");
                   PostData(1, fListHistCascade);
                   PostData(2, fCFContCascadePIDAsXiMinus);
                   PostData(3, fCFContCascadePIDAsXiPlus);
                   PostData(4, fCFContCascadePIDAsOmegaMinus);
                   PostData(5, fCFContCascadePIDAsOmegaPlus);
                   PostData(6, fCFContAsCascadeCuts);
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
                   PostData(2, fCFContCascadePIDAsXiMinus);
                   PostData(3, fCFContCascadePIDAsXiPlus);
                   PostData(4, fCFContCascadePIDAsOmegaMinus);
                   PostData(5, fCFContCascadePIDAsOmegaPlus);
                   PostData(6, fCFContAsCascadeCuts);
                   return;
              }
              if (fHasVertex == kFALSE) { //Is First event in chunk rejection: Still present!  //FIXME
                   AliWarning("Pb / This is the first event in the chunk! ... return !");
                   PostData(1, fListHistCascade);
                   PostData(2, fCFContCascadePIDAsXiMinus);
                   PostData(3, fCFContCascadePIDAsXiPlus);
                   PostData(4, fCFContCascadePIDAsOmegaMinus);
                   PostData(5, fCFContCascadePIDAsOmegaPlus);
                   PostData(6, fCFContAsCascadeCuts);
                   return;
              } 
          }
   }
   fHistCascadeMultiplicityAfterVerticesProximity->Fill(ncascades);

   // - (G) Efficiency denominator after proximity check for the vertices
   Int_t lNbMCPrimary_G = 0;   lNbMCPrimary_G = lMCstack->GetNprimary();
   Int_t lPrimaryTrackMultiplicity_G = nTrackMultiplicity;
   Double_t partEnergy_G, partPz_G, partPt_G, lRapXiMC_G;
   Int_t PDGcode_G;
   TParticle* lCurrentParticlePrimary_G = 0x0;
   AliAODMCParticle *lCurrentParticleaod_G = 0x0;
   for (Int_t iCurrentLabelStack_G = 0; iCurrentLabelStack_G < lNbMCPrimary_G; iCurrentLabelStack_G++) {
         partEnergy_G = 0.;  partPz_G = 0.;  partPt_G = 0.;  PDGcode_G = 0;
         if (fAnalysisType == "ESD") {
              lCurrentParticlePrimary_G = 0x0;  lCurrentParticlePrimary_G = lMCstack->Particle(iCurrentLabelStack_G);
              if (!lCurrentParticlePrimary_G) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack_G)) continue;
              partEnergy_G = lCurrentParticlePrimary_G->Energy();  partPz_G = lCurrentParticlePrimary_G->Pz();  partPt_G = lCurrentParticlePrimary_G->Pt();  PDGcode_G = lCurrentParticlePrimary_G->GetPdgCode();
         } else if (fAnalysisType == "AOD") {
              lCurrentParticleaod_G = 0x0;  lCurrentParticleaod_G = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack_G);
              if (!lCurrentParticleaod_G) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
              if (!lCurrentParticleaod_G->IsPhysicalPrimary()) continue;
              partEnergy_G = lCurrentParticleaod_G->E();  partPz_G = lCurrentParticleaod_G->Pz();  partPt_G = lCurrentParticleaod_G->Pt();  PDGcode_G = lCurrentParticleaod_G->GetPdgCode();
         }
         lRapXiMC_G = 0.5*TMath::Log((partEnergy_G + partPz_G) / (partEnergy_G - partPz_G + 1.e-13));
         if (PDGcode_G ==  3312)  f3dHistGenPtVsGenYvsNtracksXiMinus_G->Fill(partPt_G, lRapXiMC_G, lPrimaryTrackMultiplicity_G);
         if (PDGcode_G == -3312)  f3dHistGenPtVsGenYvsNtracksXiPlus_G->Fill(partPt_G, lRapXiMC_G, lPrimaryTrackMultiplicity_G);
         if (PDGcode_G ==  3334)  f3dHistGenPtVsGenYvsNtracksOmegaMinus_G->Fill(partPt_G, lRapXiMC_G, lPrimaryTrackMultiplicity_G);
         if (PDGcode_G == -3334)  f3dHistGenPtVsGenYvsNtracksOmegaPlus_G->Fill(partPt_G, lRapXiMC_G, lPrimaryTrackMultiplicity_G);
   }


   //=========================================================
   // 3.4) Vertex Z position selection (+ magnetic field info)
   //=========================================================
   // - Vertex coordinates: get the best primary vertex available for the event
   Double_t lBestPrimaryVtxPos[3]  = {-100.0, -100.0, -100.0};
   Double_t tPrimaryVtxPosition[3] = {-100.0, -100.0, -100.0};
   if (fAnalysisType == "ESD") {
      const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex();
      if (!lPrimaryBestESDVtx){
          AliWarning("No prim. vertex in ESD... return!");
          PostData(1, fListHistCascade);
          PostData(2, fCFContCascadePIDAsXiMinus);
          PostData(3, fCFContCascadePIDAsXiPlus);
          PostData(4, fCFContCascadePIDAsOmegaMinus);
          PostData(5, fCFContCascadePIDAsOmegaPlus);
          PostData(6, fCFContAsCascadeCuts);
          return;
      }
      lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
      // - Fill the vertex plots before any event selection on vertex position
      const AliVVertex *primaryVtx = lESDevent->GetPrimaryVertex();
      tPrimaryVtxPosition[0] = primaryVtx->GetX();
      tPrimaryVtxPosition[1] = primaryVtx->GetY();
      tPrimaryVtxPosition[2] = primaryVtx->GetZ();
      fHistPVx->Fill( tPrimaryVtxPosition[0] );
      fHistPVy->Fill( tPrimaryVtxPosition[1] );
      fHistPVz->Fill( tPrimaryVtxPosition[2] );
   } else if (fAnalysisType == "AOD") {
      const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();
      if (!lPrimaryBestAODVtx){
          AliWarning("No prim. vertex in AOD... return!");
          PostData(1, fListHistCascade);
          PostData(2, fCFContCascadePIDAsXiMinus);
          PostData(3, fCFContCascadePIDAsXiPlus);
          PostData(4, fCFContCascadePIDAsOmegaMinus);
          PostData(5, fCFContCascadePIDAsOmegaPlus);
          PostData(6, fCFContAsCascadeCuts);
          return;
      }
      lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );
      // - Fill the vertex plots before any event selection on vertex position
      const AliVVertex *primaryVtx = lAODevent->GetPrimaryVertex();
      tPrimaryVtxPosition[0] = primaryVtx->GetX();
      tPrimaryVtxPosition[1] = primaryVtx->GetY();
      tPrimaryVtxPosition[2] = primaryVtx->GetZ();
      fHistPVx->Fill( tPrimaryVtxPosition[0] );
      fHistPVy->Fill( tPrimaryVtxPosition[1] );
      fHistPVz->Fill( tPrimaryVtxPosition[2] );
   }
   // - Selection on the primary vertex Z position  
   if (fApplyEvSelZprimVtxPos) {
      if (fAnalysisType == "ESD") {
          if (TMath::Abs(lBestPrimaryVtxPos[2]) < fVtxRangeMin || TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRangeMax) {
               AliWarning("Pb / | Z position of Best Prim Vtx | out of range [Rmin, Rmax]... return !");
               PostData(1, fListHistCascade);
               PostData(2, fCFContCascadePIDAsXiMinus);
               PostData(3, fCFContCascadePIDAsXiPlus);
               PostData(4, fCFContCascadePIDAsOmegaMinus);
               PostData(5, fCFContCascadePIDAsOmegaPlus);
               PostData(6, fCFContAsCascadeCuts);
               return;
          }
      } else if (fAnalysisType == "AOD") {
          if (TMath::Abs(lBestPrimaryVtxPos[2]) < fVtxRangeMin || TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRangeMax) {
              AliWarning("Pb / | Z position of Best Prim Vtx | out of range [Rmin, Rmax]... return !");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDAsXiMinus);
              PostData(3, fCFContCascadePIDAsXiPlus);
              PostData(4, fCFContCascadePIDAsOmegaMinus);
              PostData(5, fCFContCascadePIDAsOmegaPlus);
              PostData(6, fCFContAsCascadeCuts);
              return;
          }
      }
   }
   fHistCascadeMultiplicityAfterZprimVtxPosSel->Fill(ncascades);
   // - Vertex position plots: after any event selections
   tPrimaryVtxPosition[0] = 0;
   tPrimaryVtxPosition[1] = 0;
   tPrimaryVtxPosition[2] = 0;
   if (fAnalysisType == "ESD" ) {
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

   //=================
   // - Magnetic field
   //=================
   Double_t lMagneticField = -10.;
   if      (fAnalysisType == "ESD") lMagneticField = lESDevent->GetMagneticField();
   else if (fAnalysisType == "AOD") lMagneticField = lAODevent->GetMagneticField();
   //if(TMath::Abs(lMagneticField ) < 10e-6) continue;    

   // - (H) Efficiency denominator after all event selections + QA plots on Generated particle distributions
   // - Initialisation of useful local variables		
   Int_t lPdgCodeCasc = 0, lPdgCodeBach = 0, lPdgCodeLambda = 0, lPdgCodeDghtMesV0 = 0, lPdgCodeDghtBarV0 = 0, ncascperev = 0, ncascperevtot =0, lPrimaryTrackMultiplicity = 0;	
   TH1F *lHistEtaGenCasc = 0, *lHistThetaGenCasc = 0, *lHistThetaLambda = 0, *lHistThetaBach = 0, *lHistThetaBarDghter = 0, *lHistThetaMesDghter = 0, 
        *lHistEtaLambda = 0, *lHistEtaBach = 0, *lHistEtaBarDghter = 0, *lHistEtaMesDghter = 0, *lHistPtBach = 0, *lHistPtBarDghter = 0, *lHistPtMesDghter = 0;
   TH3D *l3dHistGenPtVsGenYvsNtracks = 0;
   TH3D *l3dHistGenPtVsGenctauvsY    = 0;
   TH2D *l2dHistGenPtVsGenYFdbl      = 0;
   Double_t partEnergy, partPz, partEta, partTheta, partP, partPt, partVx, partVy, partVz, bacVx, bacVy, bacVz, partMass, lRapXiMC, lctau, lRadToDeg;
   Float_t lambdaTheta, bacTheta, dghtBarV0Theta, dghtMesV0Theta,lambdaEta, bacEta, dghtBarV0Eta, dghtMesV0Eta, bacPt, dghtBarV0Pt, dghtMesV0Pt;
   TParticle *lCurrentParticle = 0x0, *xiMC = 0x0, *mcBach = 0x0, *lDght0ofXi = 0x0, *lDght1ofXi = 0x0, *lLambda = 0x0, 
             *lBach = 0x0, *lDghtBarV0 = 0x0, *lDghtMesV0 = 0x0, *lDght0ofLambda = 0x0, *lDght1ofLambda = 0x0;
   AliAODMCParticle *lCurrentParticleaod = 0x0, *xiMCaod = 0x0, *mcBachaod = 0x0, *lDght0ofXiaod = 0x0, *lDght1ofXiaod = 0x0, *lLambdaaod = 0x0, 
                    *lBachaod = 0x0, *lDghtBarV0aod = 0x0, *lDghtMesV0aod = 0x0, *lDght0ofLambdaaod = 0x0, *lDght1ofLambdaaod = 0x0;
   Int_t lNbMCPrimary_H = 0;  lNbMCPrimary_H = lMCstack->GetNprimary();
   // - Start loop over different kind of cascades (Xi-+, Omega-+)
   for (Int_t iCascType = 1; iCascType < 5; iCascType++) { 
         ncascperev = 0;
         ncascperevtot = 0;
         lPrimaryTrackMultiplicity = nTrackMultiplicity;
         switch (iCascType) {
           case 1: // Xi-
               lPdgCodeCasc                = 3312;                                 // Xi- pdg code
               lPdgCodeBach                = -211;                                 // Pi- pdg code
               lPdgCodeLambda              = 3122;                                 // Lambda0 pdg code
               lPdgCodeDghtMesV0           = -211;                                 // Pi- pdg code
               lPdgCodeDghtBarV0           = 2212;                                 // Proton pdg code
               lHistEtaGenCasc             = fHistEtaGenCascXiMinus;               // this plot for any Xi- 
  	       lHistThetaGenCasc           = fHistThetaGenCascXiMinus;             // cascades generated within acceptance (cut in pt + theta)
               l3dHistGenPtVsGenYvsNtracks = f3dHistGenPtVsGenYvsNtracksXiMinus_H; //
               l3dHistGenPtVsGenctauvsY    = f3dHistGenPtVsGenctauvsYXiMinus_H;    //
	       l2dHistGenPtVsGenYFdbl      = f2dHistGenPtVsGenYFdblXiMinus;        //
	       lHistThetaLambda            = fHistThetaLambdaXiMinus;              // 
	       lHistThetaBach              = fHistThetaBachXiMinus;                // 
	       lHistThetaBarDghter         = fHistThetaBarDghterXiMinus;           // 
	       lHistThetaMesDghter         = fHistThetaMesDghterXiMinus;           // 
               lHistEtaLambda              = fHistEtaLambdaXiMinus;                // 
               lHistEtaBach                = fHistEtaBachXiMinus;                  // 
               lHistEtaBarDghter           = fHistEtaBarDghterXiMinus;             // 
               lHistEtaMesDghter           = fHistEtaMesDghterXiMinus;             //
	       lHistPtBach	           = fHistPtBachXiMinus;                   // 
	       lHistPtBarDghter            = fHistPtBarDghterXiMinus;              // 
	       lHistPtMesDghter            = fHistPtMesDghterXiMinus;              // 
               break; 
           case 2: // Xi+
               lPdgCodeCasc                = -3312;                               // Xi+ pdg code
               lPdgCodeBach                = 211;                                 // Pi+ pdg code
               lPdgCodeLambda              = -3122;                               // AntiLambda0 pdg code
               lPdgCodeDghtMesV0           = 211;                                 // Pi+ pdg code
               lPdgCodeDghtBarV0           = -2212;                               // AntiProton pdg code 
      	       lHistEtaGenCasc             = fHistEtaGenCascXiPlus;               // this plot for any Xi+
	       lHistThetaGenCasc           = fHistThetaGenCascXiPlus;             // cascades generated within acceptance (cut in pt + theta)
               l3dHistGenPtVsGenYvsNtracks = f3dHistGenPtVsGenYvsNtracksXiPlus_H; //
               l3dHistGenPtVsGenctauvsY    = f3dHistGenPtVsGenctauvsYXiPlus_H;    //
	       l2dHistGenPtVsGenYFdbl      = f2dHistGenPtVsGenYFdblXiPlus;        //
	       lHistThetaLambda            = fHistThetaLambdaXiPlus;              //
	       lHistThetaBach              = fHistThetaBachXiPlus;                //
	       lHistThetaBarDghter         = fHistThetaBarDghterXiPlus;           // 
	       lHistThetaMesDghter         = fHistThetaMesDghterXiPlus;           //
               lHistEtaLambda              = fHistEtaLambdaXiPlus;                //
               lHistEtaBach                = fHistEtaBachXiPlus;                  //
               lHistEtaBarDghter           = fHistEtaBarDghterXiPlus;             // 
               lHistEtaMesDghter           = fHistEtaMesDghterXiPlus;             // 
	       lHistPtBach	           = fHistPtBachXiPlus;                   //
	       lHistPtBarDghter            = fHistPtBarDghterXiPlus;              // 
	       lHistPtMesDghter            = fHistPtMesDghterXiPlus;              // 
    	       break;
           case 3: // Omega-
    	       lPdgCodeCasc                =   3334;                                  // Omega- pdg code
               lPdgCodeBach                =   -321;                                  // K- pdg code
               lPdgCodeLambda              =   3122;                                  // Lambda0 pdg code
               lPdgCodeDghtMesV0           =   -211;                                  // Pi- pdg code
               lPdgCodeDghtBarV0           =   2212;                                  // Proton pdg code
	       lHistEtaGenCasc             = fHistEtaGenCascOmegaMinus;               // this plot for any Omega+	 	
	       lHistThetaGenCasc           = fHistThetaGenCascOmegaMinus;             // cascades generated within acceptance (cut in pt + theta)
	       l2dHistGenPtVsGenYFdbl      = f2dHistGenPtVsGenYFdblOmegaMinus;        //
               l3dHistGenPtVsGenYvsNtracks = f3dHistGenPtVsGenYvsNtracksOmegaMinus_H; //
               l3dHistGenPtVsGenctauvsY    = f3dHistGenPtVsGenctauvsYOmegaMinus_H;    //
	       lHistThetaLambda            = fHistThetaLambdaOmegaMinus;              //
	       lHistThetaBach              = fHistThetaBachOmegaMinus;                //
	       lHistThetaBarDghter         = fHistThetaBarDghterOmegaMinus;           //
	       lHistThetaMesDghter         = fHistThetaMesDghterOmegaMinus;           //
               lHistEtaLambda              = fHistEtaLambdaOmegaMinus;                //
               lHistEtaBach                = fHistEtaBachOmegaMinus;                  //
               lHistEtaBarDghter           = fHistEtaBarDghterOmegaMinus;             //
               lHistEtaMesDghter           = fHistEtaMesDghterOmegaMinus;             //
	       lHistPtBach	           = fHistPtBachOmegaMinus;                   //
	       lHistPtBarDghter            = fHistPtBarDghterOmegaMinus;              //
	       lHistPtMesDghter            = fHistPtMesDghterOmegaMinus;              //
               break;
           case 4:  // Omega+
               lPdgCodeCasc                = -3334;                                  // Omega+ pdg code
               lPdgCodeBach                = 321;                                    // K+ pdg code
               lPdgCodeLambda              = -3122;                                  // AntiLambda0 pdg code
               lPdgCodeDghtMesV0           = 211;                                    // Pi+ pdg code
               lPdgCodeDghtBarV0           = -2212;                                  // AntiProton pdg code
	       lHistEtaGenCasc             = fHistEtaGenCascOmegaPlus;               // this plot for any Omega-
	       lHistThetaGenCasc           = fHistThetaGenCascOmegaPlus;             // cascades generated within acceptance (cut in pt + theta)
	       l2dHistGenPtVsGenYFdbl      = f2dHistGenPtVsGenYFdblOmegaPlus;        //
               l3dHistGenPtVsGenYvsNtracks = f3dHistGenPtVsGenYvsNtracksOmegaPlus_H; //
               l3dHistGenPtVsGenctauvsY    = f3dHistGenPtVsGenctauvsYOmegaPlus_H;    //
	       lHistThetaLambda            = fHistThetaLambdaOmegaPlus;              //
	       lHistThetaBach              = fHistThetaBachOmegaPlus;                //
	       lHistThetaBarDghter         = fHistThetaBarDghterOmegaPlus;           //
	       lHistThetaMesDghter         = fHistThetaMesDghterOmegaPlus;           //
               lHistEtaLambda              = fHistEtaLambdaOmegaPlus;                //
               lHistEtaBach                = fHistEtaBachOmegaPlus;                  //
               lHistEtaBarDghter           = fHistEtaBarDghterOmegaPlus;             //
               lHistEtaMesDghter           = fHistEtaMesDghterOmegaPlus;             //
	       lHistPtBach	           = fHistPtBachOmegaPlus;                   //
	       lHistPtBarDghter            = fHistPtBarDghterOmegaPlus;              //
	       lHistPtMesDghter            = fHistPtMesDghterOmegaPlus;              //
               break;
         }
         // - Start loop on primaries
         for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < lNbMCPrimary_H; iCurrentLabelStack++) {
               partEnergy = 0.;   partPz = 0.;  partEta = 0.;  partTheta = 0.;  partP = 0.;  partPt = 0.;
               partVx     = 0.;   partVy = 0.;  partVz = 0.;   bacVx = 0.;      bacVy = 0.;  bacVz  = 0.; partMass   = 0.;
               if ( fAnalysisType == "ESD" ) {      
                    lCurrentParticle = 0x0;  lCurrentParticle = lMCstack->Particle( iCurrentLabelStack );
                    if (!lCurrentParticle) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
                    if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) continue; 
                    if (lCurrentParticle->GetPdgCode() == lPdgCodeCasc) {	
                         partEnergy = lCurrentParticle->Energy(); partPz   = lCurrentParticle->Pz();       partEta = lCurrentParticle->Eta(); partP  = lCurrentParticle->P();  partPt = lCurrentParticle->Pt();
                         partTheta  = lCurrentParticle->Theta();  partMass = lCurrentParticle->GetMass();  partVx = lCurrentParticle->Vx();   partVy = lCurrentParticle->Vy(); partVz = lCurrentParticle->Vz();
                         if (lCurrentParticle->GetDaughter(0) >= 0) {
                              mcBach = 0x0;  mcBach = lMCstack->Particle(lCurrentParticle->GetDaughter(0));
                              if (mcBach) { bacVx = mcBach->Vx();  bacVy = mcBach->Vy();  bacVz = mcBach->Vz(); }
                         }
                    } else continue;
               } else if ( fAnalysisType == "AOD" ) {
                    lCurrentParticleaod = 0x0;  lCurrentParticleaod = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack);
                    if (!lCurrentParticleaod) { AliWarning("Generated cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n");  continue; }
                    if (!lCurrentParticleaod->IsPhysicalPrimary()) continue;  
                    if (!(lCurrentParticleaod->PdgCode() == lPdgCodeCasc)) continue;
                    partEnergy = lCurrentParticleaod->E();     partPz = lCurrentParticleaod->Pz();    partEta = lCurrentParticleaod->Eta();   partP = lCurrentParticleaod->P();    partPt = lCurrentParticleaod->Pt();   
                    partTheta = lCurrentParticleaod->Theta();  partMass = lCurrentParticleaod->M();   partVx = lCurrentParticleaod->Xv();     partVy = lCurrentParticleaod->Yv();  partVz = lCurrentParticleaod->Zv();
                    if (lCurrentParticleaod->GetDaughter(0) >= 0) {
                         mcBachaod = 0x0;  mcBachaod = (AliAODMCParticle*) arrayMC->At(lCurrentParticleaod->GetDaughter(0));
                         if (mcBachaod) {  bacVx = mcBachaod->Xv();  bacVy = mcBachaod->Yv();  bacVz = mcBachaod->Zv();  } 
                    }
               }
               ncascperevtot++;	
               lRapXiMC = 0.5*TMath::Log((partEnergy + partPz) / (partEnergy - partPz +1.e-13));
               lctau = TMath::Sqrt((partVx-bacVx)*(partVx-bacVx)+(partVy-bacVy)*(partVy-bacVy)+(partVz-bacVz)*(partVz-bacVz));
               if (partP!=0.) lctau = lctau*partMass/partP;
               else           lctau = -1.;
               lRadToDeg = 180.0/TMath::Pi();
               // - Fill histos for generated cascade not necessarily within the acceptance		
               lHistEtaGenCasc->Fill(partEta);	 
               l3dHistGenPtVsGenYvsNtracks->Fill(partPt, lRapXiMC, lPrimaryTrackMultiplicity);
               l3dHistGenPtVsGenctauvsY->Fill(partPt, lctau, lRapXiMC);
               lHistThetaGenCasc->Fill(lRadToDeg * partTheta);
               // - Check the emission of particle stays within the acceptance of the detector (cut in theta)
               if (fApplyAccCut) {if( partTheta < TMath::Pi()/4.0 || partTheta > 3.0*TMath::Pi()/4.0 ) continue;}	
               lambdaTheta = 0.;  bacTheta = 0.;  dghtBarV0Theta = 0.;  dghtMesV0Theta = 0.; lambdaEta = 0.;  bacEta = 0.;  dghtBarV0Eta = 0.;  dghtMesV0Eta = 0.;  bacPt = 0.;  dghtBarV0Pt = 0.;  dghtMesV0Pt = 0.;
               if (fAnalysisType == "ESD") { 
                    // - Cascade level
                    xiMC = 0;  xiMC = lMCstack->Particle( iCurrentLabelStack );
                    if ( xiMC->GetNDaughters() != 2) continue;
                    if ( xiMC->GetDaughter(0) < 0 )  continue;
                    if ( xiMC->GetDaughter(1) < 0 )  continue;	
                    lDght0ofXi = 0;  lDght0ofXi = lMCstack->Particle(xiMC->GetDaughter(0));
                    lDght1ofXi = 0;  lDght1ofXi = lMCstack->Particle(xiMC->GetDaughter(1));
                    lLambda = 0;  lBach = 0;
                    // -- Case 1
                    if ( lDght0ofXi->GetPdgCode() == lPdgCodeLambda && lDght1ofXi->GetPdgCode() == lPdgCodeBach ){  			
                           lLambda = lDght0ofXi;   // dghter0 = Lambda
                           lBach   = lDght1ofXi;   // dghter1 = Pi-
                    }  		
                    // -- Case 2
                    else if ( lDght0ofXi->GetPdgCode() == lPdgCodeBach && lDght1ofXi->GetPdgCode() == lPdgCodeLambda ){ 
                           lBach   = lDght0ofXi;   // dghter0 = Pi-
	                   lLambda = lDght1ofXi;   // dghter1 = Lambda
                    }
	            // -- Case 3	
                    else continue;
                    // - Check the emission of particle stays within the acceptance of the detector (cut in pt + theta)
                    if (fApplyAccCut) { 
                         if( lLambda->Theta() < TMath::Pi()/4.0 || lLambda->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
                         if( lBach->Theta() < TMath::Pi()/4.0   || lBach->Theta() > 3.0*TMath::Pi()/4.0 )   continue;
                         if( lBach->Pt() < 0.150 ) continue; //FIXME: maybe tuned for Xi but not for K- from Omega ...
                    } 				
                    // - V0 level
                    lDghtBarV0 = 0;  lDghtMesV0 = 0;
                    if( lLambda->GetNDaughters() != 2 )  continue;
                    if( lLambda->GetDaughter(0) < 0 )    continue;
                    if( lLambda->GetDaughter(1) < 0 )    continue;
                    lDght0ofLambda = 0;  lDght0ofLambda = lMCstack->Particle(lLambda->GetDaughter(0));
                    lDght1ofLambda = 0;  lDght1ofLambda = lMCstack->Particle(lLambda->GetDaughter(1));
                    // -- Case 1
                    if ( lDght0ofLambda->GetPdgCode() == lPdgCodeDghtBarV0 && lDght1ofLambda->GetPdgCode() == lPdgCodeDghtMesV0 ) { 
                           lDghtBarV0 = lDght0ofLambda;   // dghter0 = Proton
                           lDghtMesV0 = lDght1ofLambda;   // dghter1 = Pi-
                    }  		
                    // -- Case 2
                    else if ( lDght0ofLambda->GetPdgCode() == lPdgCodeDghtMesV0 && lDght1ofLambda->GetPdgCode() == lPdgCodeDghtBarV0 ) { 
                           lDghtMesV0 = lDght0ofLambda;   // dghter0 = Pi-
	                   lDghtBarV0 = lDght1ofLambda;   // dghter1 = Proton
                    }		
                    // -- Case 3
                    else continue;
                    // - Check the emission of particle stays within the acceptance of the detector
                    if (fApplyAccCut) { 
                         if( lDghtBarV0->Theta() < TMath::Pi()/4.0  ||  lDghtBarV0->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
                         if( lDghtMesV0->Theta() < TMath::Pi()/4.0  ||  lDghtMesV0->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
                         if( lDghtBarV0->Pt() < 0.250 ) continue;
                         if( lDghtMesV0->Pt() < 0.150 ) continue;
                    }
                    lambdaTheta    = lLambda->Theta();
                    bacTheta       = lBach->Theta();
                    dghtBarV0Theta = lDghtBarV0->Theta(); 			
                    dghtMesV0Theta = lDghtMesV0->Theta();
                    lambdaEta      = lLambda->Eta();
                    bacEta         = lBach->Eta();
                    dghtBarV0Eta   = lDghtBarV0->Eta();
                    dghtMesV0Eta   = lDghtMesV0->Eta();
                    bacPt          = lBach->Pt();
                    dghtBarV0Pt    = lDghtBarV0->Pt();
                    dghtMesV0Pt    = lDghtMesV0->Pt();
               } else if ( fAnalysisType == "AOD") {
                    // - Cascade level
                    xiMCaod = 0x0;  xiMCaod = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack);
                    if (xiMCaod->GetNDaughters() != 2) continue;
                    if (xiMCaod->GetDaughter(0) < 0 )  continue;
                    if (xiMCaod->GetDaughter(1) < 0 )  continue;
                    lDght0ofXiaod = 0x0;  lDght0ofXiaod = (AliAODMCParticle*) arrayMC->At(xiMCaod->GetDaughter(0));
                    lDght1ofXiaod = 0x0;  lDght1ofXiaod = (AliAODMCParticle*) arrayMC->At(xiMCaod->GetDaughter(1));
                    lLambdaaod = 0x0;   lBachaod = 0x0;
                    // -- Case 1
                    if ( lDght0ofXiaod->PdgCode() == lPdgCodeLambda  &&  lDght1ofXiaod->PdgCode() == lPdgCodeBach ){ 
                            lLambdaaod = lDght0ofXiaod;   // dghter0 = Lambda
                            lBachaod   = lDght1ofXiaod;   // dghter1 = Pi-
                    }
                    // -- Case 2
                    else if ( lDght0ofXiaod->PdgCode() == lPdgCodeBach && lDght1ofXiaod->PdgCode() == lPdgCodeLambda ){
                            lBachaod   = lDght0ofXiaod;   // dghter0 = Pi
                            lLambdaaod = lDght1ofXiaod;   //dghter1 = Lambda
                    }
                    // -- Case 3
                    else continue;
                    // - Check the emission of particle stays within the acceptance of the detector (cut in pt + theta)
                    if (fApplyAccCut) {
                          if( lLambdaaod->Theta() < TMath::Pi()/4.0  ||    lLambdaaod->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
                          if( lBachaod->Theta() < TMath::Pi()/4.0    ||    lBachaod->Theta() > 3.0*TMath::Pi()/4.0 )   continue;
                          if( lBachaod->Pt() < 0.150 ) continue; //FIXME : maybe tuned for Xi but not for K- from Omega ...
                    }
                    // - V0 level 
                    lDghtBarV0aod = 0x0;  lDghtMesV0aod = 0x0;
                    if( lLambdaaod->GetNDaughters() != 2 )  continue;
                    if( lLambdaaod->GetDaughter(0) < 0 )    continue;
                    if( lLambdaaod->GetDaughter(1) < 0 )    continue;
                    lDght0ofLambdaaod = 0x0;  lDght0ofLambdaaod = (AliAODMCParticle*) arrayMC->At(lLambdaaod->GetDaughter(0));
                    lDght1ofLambdaaod = 0x0;  lDght1ofLambdaaod = (AliAODMCParticle*) arrayMC->At(lLambdaaod->GetDaughter(1));
                    // -- Case 1
                    if ( lDght0ofLambdaaod->PdgCode() == lPdgCodeDghtBarV0 && lDght1ofLambdaaod->PdgCode() == lPdgCodeDghtMesV0 ) { 
                            lDghtBarV0aod = lDght0ofLambdaaod;   // dghter0 = Proton
                            lDghtMesV0aod = lDght1ofLambdaaod;   // dghter1 = Pi-
                    } 
                    // -- Case 2
                    else if ( lDght0ofLambdaaod->PdgCode() == lPdgCodeDghtMesV0 && lDght1ofLambdaaod->PdgCode() == lPdgCodeDghtBarV0 ) { 
                            lDghtMesV0aod = lDght0ofLambdaaod;   // dghter0 = Pi-
                            lDghtBarV0aod = lDght1ofLambdaaod;   // dghter1 = proton
                    } 
                    // -- Case 3
                    else continue;
                    // - Check the emission of particle stays within the acceptance of the detector
                    if (fApplyAccCut) {
                         if( lDghtBarV0aod->Theta() < TMath::Pi()/4.0  ||  lDghtBarV0aod->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
                         if( lDghtMesV0aod->Theta() < TMath::Pi()/4.0  ||  lDghtMesV0aod->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
                         if( lDghtBarV0aod->Pt() < 0.250 ) continue;
                         if( lDghtMesV0aod->Pt() < 0.150 ) continue;
                    }
                    lambdaTheta    = lLambdaaod->Theta();
                    bacTheta       = lBachaod->Theta();
                    dghtBarV0Theta = lDghtBarV0aod->Theta();
                    dghtMesV0Theta = lDghtMesV0aod->Theta();
                    lambdaEta      = lLambdaaod->Eta();
                    bacEta         = lBachaod->Eta();
                    dghtBarV0Eta   = lDghtBarV0aod->Eta();
                    dghtMesV0Eta   = lDghtMesV0aod->Eta();
                    bacPt          = lBachaod->Pt();
                    dghtBarV0Pt    = lDghtBarV0aod->Pt();
                    dghtMesV0Pt    = lDghtMesV0aod->Pt();
               }
               // - Filling histos for findable cascades (findable = cascade within the acceptance)
               lHistThetaLambda->Fill( lRadToDeg * lambdaTheta );
               lHistThetaBach->Fill( lRadToDeg * bacTheta );
               lHistThetaBarDghter->Fill( lRadToDeg * dghtBarV0Theta );
               lHistThetaMesDghter->Fill( lRadToDeg * dghtMesV0Theta );
               if (partEta >= 1.2) lHistEtaLambda->Fill( lambdaEta );
               if (partEta >= 1.2) lHistEtaBach->Fill( bacEta );
               if (partEta >= 1.2) lHistEtaBarDghter->Fill( dghtBarV0Eta );
               if (partEta >= 1.2) lHistEtaMesDghter->Fill( dghtMesV0Eta );
               lHistPtBach->Fill( bacPt );
               lHistPtBarDghter->Fill( dghtBarV0Pt );
               lHistPtMesDghter->Fill( dghtMesV0Pt );
               l2dHistGenPtVsGenYFdbl->Fill( partPt, lRapXiMC );
               // - Increase the counter of cascades
               ncascperev++;
         } // This is the end of the loop on primaries
         // - Fill the histo containing the number of cascades per event before and after the fApplyAccCut
         if (iCascType == 1) {
              fHistnXiMinusPerEv->Fill(ncascperev);
              fHistnXiMinusPerEvTot->Fill(ncascperevtot);
         }
         if (iCascType == 2) {
              fHistnXiPlusPerEv->Fill(ncascperev);
              fHistnXiPlusPerEvTot->Fill(ncascperevtot);
         }
         if (iCascType == 3) {
              fHistnOmegaMinusPerEv->Fill(ncascperev);
              fHistnOmegaMinusPerEvTot->Fill(ncascperevtot);
         }
         if (iCascType == 4) {
              fHistnOmegaPlusPerEv->Fill(ncascperev);
              fHistnOmegaPlusPerEvTot->Fill(ncascperevtot);
         }
         // - Re-initialisation of the local THF pointers
         lHistEtaGenCasc = 0x0;      lHistThetaGenCasc = 0x0;    l2dHistGenPtVsGenYFdbl = 0x0;    lHistThetaLambda = 0x0;   lHistThetaBach = 0x0;   
         lHistThetaBarDghter = 0x0;  lHistThetaMesDghter = 0x0;  lHistPtBach = 0x0;               lHistPtBarDghter = 0x0;   lHistPtMesDghter = 0x0;	
   } // end of loop over the different types of cascades (Xi-+, Omega-+)
 	
 
 
   //-----------------------------------------	
   // - Loop over the reconstructed candidates
   //-----------------------------------------
   Int_t nAssoXiMinus     = 0;
   Int_t nAssoXiPlus      = 0;
   Int_t nAssoOmegaMinus  = 0;
   Int_t nAssoOmegaPlus   = 0;
   Int_t lPosTPCClusters  = 0;
   Int_t lNegTPCClusters  = 0;
   Int_t lBachTPCClusters = 0;
   Int_t lblBachForPID; 
   Int_t lblPosV0Dghter; 
   Int_t lblNegV0Dghter; 
   Int_t lblMotherPosV0Dghter;
   Int_t lblMotherNegV0Dghter; 
   Int_t lblBach;
   Int_t lblGdMotherPosV0Dghter; 
   Int_t lblGdMotherNegV0Dghter; 
   Int_t lblMotherBach;
   Double_t lPosTPCFindClusters         = -1;
   Double_t lNegTPCFindClusters         = -1;
   Double_t lBachTPCFindClusters        = -1;
   Double_t lDcaXiDaughters = 0., lDcaBachToPrimVertexXi = 0., lXiCosineOfPointingAngle = 0., lXiRadius = 0., lInvMassLambdaAsCascDghter = 0., lDcaV0DaughtersXi = 0., lV0CosineOfPointingAngleXi = 0., 
            lV0CosineOfPointingAngle = 0., lV0RadiusXi = 0., lDcaV0ToPrimVertexXi = 0., lDcaPosToPrimVertexXi = 0., lDcaNegToPrimVertexXi = 0., lChargeXi = 0. , lV0mom = 0., lmcPt = 0., lmcRapCasc = 0., 
            lmcEta = 0., lmcTransvRadius = 0., lrecoPt = 0., lrecoTransvRadius = 0., lDeltaPhiMcReco = 0., lBachTransvMom = 0., lpTrackTransvMom = 0., lnTrackTransvMom = 0., lmcPtPosV0Dghter = 0., 
            lmcPtNegV0Dghter = 0., lrecoP = 0., lmcPtBach = 0., cascadeMass = 0., ppionBach, pkaonBach, lInvMassXiMinus, lInvMassXiPlus, lInvMassOmegaMinus, lInvMassOmegaPlus, lV0quality, pproton, ppion;
   Bool_t   lIsPosInXiProton, lIsPosInXiPion, lIsPosInOmegaProton, lIsPosInOmegaPion, lIsNegInXiProton, lIsNegInXiPion, lIsNegInOmegaProton, lIsNegInOmegaPion, lIsBachelorKaon, lIsBachelorPion, 
            lIsBachelorKaonForTPC, lIsBachelorPionForTPC, lIsNegPionForTPC, lIsPosPionForTPC, lIsNegProtonForTPC, lIsPosProtonForTPC, lIsBachelorMCPiMinus, lIsBachelorMCPiPlus, lIsBachelorMCKMinus, 
            lIsBachelorMCKPlus, lAssoXiMinus, lAssoXiPlus, lAssoOmegaMinus, lAssoOmegaPlus;
   Float_t  etaBach, etaPos, etaNeg, decayCascX, decayCascY, distV0Xi, lctauV0, distTV0Xi;
   AliESDcascade *xiESD;
   UInt_t lIdxPosXi, lIdxNegXi, lBachIdx;
   AliESDtrack *pTrackXi, *nTrackXi, *bachTrackXi;
   ULong_t pStatus, nStatus, bachStatus;
   TParticle *mcBachForPID, *mcPosV0Dghter, *mcNegV0Dghter, *mcMotherPosV0Dghter, *mcMotherNegV0Dghter, *mcGdMotherPosV0Dghter, *mcGdMotherNegV0Dghter, *mcMotherBach;
   const AliAODcascade *xiAOD;
   AliAODTrack *pTrackXiaod, *nTrackXiaod, *bachTrackXiaod;
   AliAODMCParticle *mcBachForPIDaod, *mcPosV0Dghteraod, *mcNegV0Dghteraod, *mcMotherPosV0Dghteraod, *mcMotherNegV0Dghteraod, *mcGdMotherPosV0Dghteraod, *mcGdMotherNegV0Dghteraod, *mcMotherBachaod;
   Double_t xiMomX, xiMomY, xiMomZ, lBachMomX, lBachMomY, lV0momX, lV0momY, lV0momZ, lV0NMomX, lV0NMomY, lV0PMomX, lV0PMomY;
   Float_t lambdaMass = 1.115683; // PDG mass 

   // - Get the number of cascades
   ncascades = 0;
   if      ( fAnalysisType == "ESD" ) { ncascades = lESDevent->GetNumberOfCascades(); }
   else if ( fAnalysisType == "AOD" ) { ncascades = lAODevent->GetNumberOfCascades(); }


   //-------------------------------
   // - Beginning of the Cascade Loop
   for (Int_t iXi = 0; iXi < ncascades; iXi++) {

        lIsPosInXiProton      = kFALSE;
        lIsPosInXiPion        = kFALSE;
        lIsPosInOmegaProton   = kFALSE;
        lIsPosInOmegaPion     = kFALSE;
        lIsNegInXiProton      = kFALSE;
        lIsNegInXiPion        = kFALSE;
        lIsNegInOmegaProton   = kFALSE;
        lIsNegInOmegaPion     = kFALSE;
        lIsBachelorKaon       = kFALSE;
        lIsBachelorPion       = kFALSE;
        lIsBachelorKaonForTPC = kFALSE;
        lIsBachelorPionForTPC = kFALSE;
        lIsNegPionForTPC      = kFALSE;
        lIsPosPionForTPC      = kFALSE;
        lIsNegProtonForTPC    = kFALSE;
        lIsPosProtonForTPC    = kFALSE;
        Double_t lPriorsGuessXi[14]    = {0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // - Combined PID - Reasonable guess for the priors for the cascade track sample (e-, mu, pi, K, p)
        Double_t lPriorsGuessOmega[14] = {0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // - Combined PID - Reasonable guess for the priors for the cascade track sample (e-, mu, pi, K, p)
        ppionBach = 0.0, pkaonBach = 0.0;
        lIsBachelorMCPiMinus  = kFALSE;
        lIsBachelorMCPiPlus   = kFALSE;
        lIsBachelorMCKMinus   = kFALSE;
        lIsBachelorMCKPlus    = kFALSE;
        lInvMassXiMinus    = 0.;
        lInvMassXiPlus     = 0.;
        lInvMassOmegaMinus = 0.;
        lInvMassOmegaPlus  = 0.;
        lAssoXiMinus    = kFALSE;
        lAssoXiPlus     = kFALSE;
        lAssoOmegaMinus = kFALSE;
        lAssoOmegaPlus  = kFALSE;
        etaBach = 0.;
        etaPos  = 0.;
        etaNeg  = 0.;
        lDcaXiDaughters          = -1;
        lDcaBachToPrimVertexXi   = -1;
        lXiCosineOfPointingAngle = -1;
        Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
        lInvMassLambdaAsCascDghter = 0.;
        lDcaV0DaughtersXi          = -1;
        lV0CosineOfPointingAngleXi = -1;
        lV0CosineOfPointingAngle   = -1;
        Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; 
        lDcaV0ToPrimVertexXi       = -1;
        lDcaPosToPrimVertexXi      = -1;
        lDcaNegToPrimVertexXi      = -1;
        lChargeXi                  = -1;
        lblBachForPID = 0;
        Double_t nV0mom[3] = {0. ,0. ,0. };
        Double_t pV0mom[3] = {0. ,0. ,0. };
        Double_t lInnerWallMomCascDghters[3] = {-100., -100., -100.};


        if ( fAnalysisType == "ESD" ) {		

             // - Load the cascade
             xiESD = 0x0;  xiESD = lESDevent->GetCascade(iXi);
	     if (!xiESD) continue;
	     // - Connection to daughter tracks of the current cascade		
             lIdxPosXi = 0x0;  lIdxPosXi = (UInt_t) TMath::Abs( xiESD->GetPindex() );
             lIdxNegXi = 0x0;  lIdxNegXi = (UInt_t) TMath::Abs( xiESD->GetNindex() );
             lBachIdx = 0x0;  lBachIdx  = (UInt_t) TMath::Abs( xiESD->GetBindex() );
             // - Rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
             if(lBachIdx == lIdxNegXi) {  AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!");   continue;  }
             if(lBachIdx == lIdxPosXi) {  AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!");   continue;  }
             // - Get the daughter tracks
	     pTrackXi = 0x0;  pTrackXi    = lESDevent->GetTrack( lIdxPosXi );
	     nTrackXi = 0x0;  nTrackXi    = lESDevent->GetTrack( lIdxNegXi );
	     bachTrackXi = 0x0;  bachTrackXi = lESDevent->GetTrack( lBachIdx  );
	     if (!pTrackXi || !nTrackXi || !bachTrackXi ) {  AliWarning("ERROR: Could not retrieve one of the 3 daughter tracks of the cascade ..."); continue;  }
             // - Rejection of a poor quality tracks
             if(fTrackQualityCutTPCrefit){  // - Poor quality related to TPCrefit
                 pStatus = 0x0;  pStatus    = pTrackXi->GetStatus();
                 nStatus = 0x0;  nStatus    = nTrackXi->GetStatus();
                 bachStatus = 0x0;  bachStatus = bachTrackXi->GetStatus();
                 if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
                 if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
                 if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
             } 
             if(fTrackQualityCutnTPCcls){ 
                lPosTPCClusters = 0;     lPosTPCClusters  = pTrackXi->GetTPCClusterInfo(2,1);      //GetTPCNcls();
                lNegTPCClusters = 0;     lNegTPCClusters  = nTrackXi->GetTPCClusterInfo(2,1);      //->GetTPCNcls();
                lBachTPCClusters = 0;    lBachTPCClusters = bachTrackXi->GetTPCClusterInfo(2,1);   //->GetTPCNcls();
                lPosTPCFindClusters   = pTrackXi->GetTPCNclsF();                  // New
                lNegTPCFindClusters   = nTrackXi->GetTPCNclsF();                  // New
                lBachTPCFindClusters  = bachTrackXi->GetTPCNclsF();               // New
                // - Poor quality related to TPC clusters
                if(lPosTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Pos. track has less than 80 TPC clusters ... continue!"); continue; }
                if(lNegTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Neg. track has less than 80 TPC clusters ... continue!"); continue; }
                if(lBachTPCClusters < fMinnTPCcls) { AliWarning("Pb / Bach.   track has less than 80 TPC clusters ... continue!"); continue; }
                // - Poor quality related to clusters/findable
                if( lPosTPCFindClusters <= 0 || lNegTPCFindClusters <= 0 || lBachTPCFindClusters ) { AliWarning("Pb / Number of findable cluster <= 0 ... continue!"); continue; }
                if ((lPosTPCClusters/lPosTPCFindClusters)    < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / V0 Pos. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
                if ((lNegTPCClusters/lNegTPCFindClusters)    < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / V0 Neg. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
                if ((lBachTPCClusters/lBachTPCFindClusters)  < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / Bach. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
             }

             const AliExternalTrackParam *pExtTrack    = pTrackXi->GetInnerParam();
             const AliExternalTrackParam *nExtTrack    = nTrackXi->GetInnerParam();
             const AliExternalTrackParam *bachExtTrack = bachTrackXi->GetInnerParam();
             if (pExtTrack)    lInnerWallMomCascDghters[0] = pExtTrack->GetP() * pExtTrack->Charge();
             if (nExtTrack)    lInnerWallMomCascDghters[1] = nExtTrack->GetP() * nExtTrack->Charge();
             if (bachExtTrack) lInnerWallMomCascDghters[2] = bachExtTrack->GetP() * bachExtTrack->Charge();


             etaPos  = pTrackXi->Eta();
             etaNeg  = nTrackXi->Eta();
             etaBach = bachTrackXi->Eta();
	     // - Info over reconstructed cascades
	     lV0quality = 0.;
	     if( bachTrackXi->Charge() < 0 ) {
                  //Calculate the effective mass of the Xi- candidate: Xi- hyp. (pdg code 3312)
	          lV0quality = 0.;
		  xiESD->ChangeMassHypothesis(lV0quality , 3312);  	
		  lInvMassXiMinus = xiESD->GetEffMassXi();
                  //Calculate the effective mass of the Xi- candidate: Omega- hyp. (pdg code 3334)
		  lV0quality = 0.;
		  xiESD->ChangeMassHypothesis(lV0quality , 3334); 	
		  lInvMassOmegaMinus = xiESD->GetEffMassXi();
                  //Back to "default" hyp. (Xi-)	
		  lV0quality = 0.;
		  xiESD->ChangeMassHypothesis(lV0quality , 3312);
	     }
	     if( bachTrackXi->Charge() >  0 ){
                  //Calculate the effective mass of the Xi- candidate: Xi+ hyp. (pdg code -3312)
	          lV0quality = 0.;
		  xiESD->ChangeMassHypothesis(lV0quality , -3312); 	
		  lInvMassXiPlus = xiESD->GetEffMassXi();
                  //Calculate the effective mass of the Xi- candidate: Omega+ hyp. (pdg code -3334)
                  lV0quality = 0.;
                  xiESD->ChangeMassHypothesis(lV0quality , -3334); 	
                  lInvMassOmegaPlus = xiESD->GetEffMassXi();
		  //Back to "default" hyp. (Xi-)
                  lV0quality = 0.;
                  xiESD->ChangeMassHypothesis(lV0quality , -3312);
	     }
             lDcaXiDaughters            = xiESD->GetDcaXiDaughters();
             lDcaBachToPrimVertexXi     = TMath::Abs( bachTrackXi->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lMagneticField) );
             lXiCosineOfPointingAngle   = xiESD->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
             xiESD->GetXYZcascade( lPosXi[0], lPosXi[1], lPosXi[2] );  
             lInvMassLambdaAsCascDghter = xiESD->GetEffMass();
             lDcaV0DaughtersXi          = xiESD->GetDcaV0Daughters();
             lV0CosineOfPointingAngleXi = xiESD->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
             lV0CosineOfPointingAngle   = xiESD->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2]); 
             xiESD->GetXYZ( lPosV0Xi[0], lPosV0Xi[1], lPosV0Xi[2] );
             lDcaV0ToPrimVertexXi       = xiESD->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
             lDcaPosToPrimVertexXi      = TMath::Abs( pTrackXi->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lMagneticField) );
             lDcaNegToPrimVertexXi      = TMath::Abs( nTrackXi->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lMagneticField) );
	     lChargeXi                  = xiESD->Charge();
	     // - PID Information
	        // - Combined V0-positive-daughter PID
	        AliPID pPidXi;     pPidXi.SetPriors(lPriorsGuessXi);
	        AliPID pPidOmega;  pPidOmega.SetPriors(lPriorsGuessOmega);		
	        if( pTrackXi->IsOn(AliESDtrack::kESDpid) ){  
	             Double_t r[10] = {0.}; 
                     pTrackXi->GetESDpid(r);
		     pPidXi.SetProbabilities(r);
		     pPidOmega.SetProbabilities(r);		
		     // Check if the V0 positive track is a proton (case for Xi-)
		     pproton = 0.;  pproton = pPidXi.GetProbability(AliPID::kProton);
		     if (pproton > pPidXi.GetProbability(AliPID::kElectron) &&  pproton > pPidXi.GetProbability(AliPID::kMuon)  &&
		         pproton > pPidXi.GetProbability(AliPID::kPion)     &&  pproton > pPidXi.GetProbability(AliPID::kKaon)  ) lIsPosInXiProton = kTRUE;
		     // Check if the V0 positive track is a pi+ (case for Xi+)
		     ppion = 0.;  ppion = pPidXi.GetProbability(AliPID::kPion);
		     if (ppion > pPidXi.GetProbability(AliPID::kElectron) &&  ppion > pPidXi.GetProbability(AliPID::kMuon) &&
		         ppion > pPidXi.GetProbability(AliPID::kKaon)     &&  ppion > pPidXi.GetProbability(AliPID::kProton)  ) lIsPosInXiPion = kTRUE;
		     // Check if the V0 positive track is a proton (case for Omega-)
		     pproton = 0.;  pproton = pPidOmega.GetProbability(AliPID::kProton);
		     if (pproton > pPidOmega.GetProbability(AliPID::kElectron) &&  pproton > pPidOmega.GetProbability(AliPID::kMuon) &&
		         pproton > pPidOmega.GetProbability(AliPID::kPion)     &&  pproton > pPidOmega.GetProbability(AliPID::kKaon)  ) lIsPosInOmegaProton = kTRUE;
	 	     // Check if the V0 positive track is a pi+ (case for Omega+)
		     ppion = 0.;  ppion = pPidOmega.GetProbability(AliPID::kPion);
		     if (ppion > pPidOmega.GetProbability(AliPID::kElectron) &&  ppion > pPidOmega.GetProbability(AliPID::kMuon) &&
		         ppion > pPidOmega.GetProbability(AliPID::kKaon)     &&  ppion > pPidOmega.GetProbability(AliPID::kProton)  ) lIsPosInOmegaPion = kTRUE;
	        }		
	        // - Combined V0-negative-daughter PID
	        AliPID nPidXi;    nPidXi.SetPriors( lPriorsGuessXi );
	        AliPID nPidOmega; nPidOmega.SetPriors( lPriorsGuessOmega );		
	        if( nTrackXi->IsOn(AliESDtrack::kESDpid) ) {  
     	             Double_t r[10] = {0.}; 
                     nTrackXi->GetESDpid(r);
	             nPidXi.SetProbabilities(r);
		     nPidOmega.SetProbabilities(r);
		     // Check if the V0 negative track is a pi- (case for Xi-)
		     ppion = 0.;  ppion = nPidXi.GetProbability(AliPID::kPion);
		     if (ppion > nPidXi.GetProbability(AliPID::kElectron) &&  ppion > nPidXi.GetProbability(AliPID::kMuon) &&
		         ppion > nPidXi.GetProbability(AliPID::kKaon)     &&  ppion > nPidXi.GetProbability(AliPID::kProton)  ) lIsNegInXiPion = kTRUE;
		     // Check if the V0 negative track is an anti-proton (case for Xi+)
		     pproton = 0.;  pproton = nPidXi.GetProbability(AliPID::kProton);
		     if (pproton > nPidXi.GetProbability(AliPID::kElectron) &&  pproton > nPidXi.GetProbability(AliPID::kMuon) &&
		         pproton > nPidXi.GetProbability(AliPID::kPion)     &&  pproton > nPidXi.GetProbability(AliPID::kKaon)  ) lIsNegInXiProton = kTRUE;
		     // Check if the V0 negative track is a pi- (case for Omega-)
		     ppion = 0.;  ppion = nPidOmega.GetProbability(AliPID::kPion);
		     if (ppion > nPidOmega.GetProbability(AliPID::kElectron) &&  ppion > nPidOmega.GetProbability(AliPID::kMuon) &&
		         ppion > nPidOmega.GetProbability(AliPID::kKaon)     &&  ppion > nPidOmega.GetProbability(AliPID::kProton)  ) lIsNegInOmegaPion = kTRUE;
		     // Check if the V0 negative track is an anti-proton (case for Omega+)
		     pproton = 0.;  pproton = nPidOmega.GetProbability(AliPID::kProton);
		     if (pproton > nPidOmega.GetProbability(AliPID::kElectron) &&  pproton > nPidOmega.GetProbability(AliPID::kMuon) &&
		         pproton > nPidOmega.GetProbability(AliPID::kPion)     &&  pproton > nPidOmega.GetProbability(AliPID::kKaon) ) lIsNegInOmegaProton = kTRUE;
	        }
	        // - Combined bachelor PID
	        AliPID bachPidXi;    bachPidXi.SetPriors( lPriorsGuessXi );
	        AliPID bachPidOmega; bachPidOmega.SetPriors( lPriorsGuessOmega );
	        if ( bachTrackXi->IsOn(AliESDtrack::kESDpid) ) {  
                     Double_t r[10] = {0.}; 
                     bachTrackXi->GetESDpid(r);
                     bachPidXi.SetProbabilities(r);
                     bachPidOmega.SetProbabilities(r);
                     // Check if the bachelor track is a pion
                     ppionBach = bachPidXi.GetProbability(AliPID::kPion);
                     if (ppionBach > bachPidXi.GetProbability(AliPID::kElectron) &&  ppionBach > bachPidXi.GetProbability(AliPID::kMuon) &&
                         ppionBach > bachPidXi.GetProbability(AliPID::kKaon)     &&  ppionBach > bachPidXi.GetProbability(AliPID::kProton) ) lIsBachelorPion = kTRUE;
                     // Check if the bachelor track is a kaon
                     pkaonBach = bachPidOmega.GetProbability(AliPID::kKaon);
                     if (pkaonBach > bachPidOmega.GetProbability(AliPID::kElectron) &&  pkaonBach > bachPidOmega.GetProbability(AliPID::kMuon) &&
                         pkaonBach > bachPidOmega.GetProbability(AliPID::kPion)     &&  pkaonBach > bachPidOmega.GetProbability(AliPID::kProton) ) lIsBachelorKaon = kTRUE;	
                }
	        // - 4-sigma bands on Bethe-Bloch curve
                // Bachelor
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < fTPCPIDsigma) lIsBachelorKaonForTPC = kTRUE;
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < fTPCPIDsigma) lIsBachelorPionForTPC = kTRUE;
                // Negative V0 daughter
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < fTPCPIDsigma) lIsNegPionForTPC   = kTRUE;
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < fTPCPIDsigma) lIsNegProtonForTPC = kTRUE;
                // Positive V0 daughter
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < fTPCPIDsigma) lIsPosPionForTPC   = kTRUE;
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < fTPCPIDsigma) lIsPosProtonForTPC = kTRUE;
	        // - PID probability vs Pt(Bach)
                lblBachForPID = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
                mcBachForPID = 0x0; mcBachForPID  = lMCstack->Particle( lblBachForPID );
                lmcPtBach = mcBachForPID->Pt();
	        // - MC perfect PID
                if (mcBachForPID->GetPdgCode() == -211) lIsBachelorMCPiMinus = kTRUE;
                if (mcBachForPID->GetPdgCode() ==  211) lIsBachelorMCPiPlus  = kTRUE;
                if (mcBachForPID->GetPdgCode() == -321) lIsBachelorMCKMinus  = kTRUE;
                if (mcBachForPID->GetPdgCode() ==  321) lIsBachelorMCKPlus   = kTRUE;
	     // -  MC association (care : lots of "continue;" below this line)
	     if (fDebug > 5) cout<< "MC EventNumber: "<<lMCevent->Header()->GetEvent()<<" / MC event Number in Run : "<<lMCevent->Header()->GetEventNrInRun()<<endl;
	     // - Level of the V0 daughters
             lblPosV0Dghter = 0;  lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );  
             lblNegV0Dghter = 0;  lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );  		
             mcPosV0Dghter = 0x0; mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
             mcNegV0Dghter = 0x0; mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
             // - Level of the cascade daughters	
             lblMotherPosV0Dghter = 0.;  lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother(); 
             lblMotherNegV0Dghter = 0.;  lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();
             if (lblMotherPosV0Dghter != lblMotherNegV0Dghter) continue; // same mother
             if (lblMotherPosV0Dghter < 0) continue;                     // this particle is primary, no mother   
             if (lblMotherNegV0Dghter < 0) continue;                     // this particle is primary, no mother
	                              	                                 // mothers = Lambda candidate ... a priori
             mcMotherPosV0Dghter = 0x0; mcMotherPosV0Dghter = lMCstack->Particle( lblMotherPosV0Dghter );
             mcMotherNegV0Dghter = 0x0; mcMotherNegV0Dghter = lMCstack->Particle( lblMotherNegV0Dghter );  // MN: redundant?? already checked that labels are the same...-->same part from stack
             lblBach = 0; lblBach = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
             mcBach = 0x0; mcBach = lMCstack->Particle( lblBach );	
             // - Level of cascade candidate
             lblGdMotherPosV0Dghter = 0; lblGdMotherPosV0Dghter = mcMotherPosV0Dghter->GetFirstMother() ;
             lblGdMotherNegV0Dghter = 0; lblGdMotherNegV0Dghter = mcMotherNegV0Dghter->GetFirstMother() ;
             if(lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter) continue;
             if(lblGdMotherPosV0Dghter < 0) continue;                  // primary lambda ...   
             if(lblGdMotherNegV0Dghter < 0) continue;                  // primary lambda ...   				
                                                                       // Gd mothers = cascade candidate ... a priori
             mcGdMotherPosV0Dghter = 0x0; mcGdMotherPosV0Dghter = lMCstack->Particle( lblGdMotherPosV0Dghter );
             mcGdMotherNegV0Dghter = 0x0; mcGdMotherNegV0Dghter = lMCstack->Particle( lblGdMotherNegV0Dghter );				
             lblMotherBach = 0; lblMotherBach = (Int_t) TMath::Abs( mcBach->GetFirstMother() );  
             if(lblMotherBach != lblGdMotherPosV0Dghter) continue; //same mother for bach and V0 daughters
             mcMotherBach = 0x0; mcMotherBach = lMCstack->Particle( lblMotherBach );
             // - Check if cascade is primary
             if (!(lMCstack->IsPhysicalPrimary(lblMotherBach))) continue;  
             // - Manage boolean for association
             if      (mcMotherBach->GetPdgCode() == 3312  && mcGdMotherPosV0Dghter->GetPdgCode() == 3312  && mcGdMotherNegV0Dghter->GetPdgCode() == 3312)  {lAssoXiMinus    = kTRUE; cascadeMass = 1.321; nAssoXiMinus++;}
             else if (mcMotherBach->GetPdgCode() == -3312 && mcGdMotherPosV0Dghter->GetPdgCode() == -3312 && mcGdMotherNegV0Dghter->GetPdgCode() == -3312) {lAssoXiPlus     = kTRUE; cascadeMass = 1.321; nAssoXiPlus++;}
             else if (mcMotherBach->GetPdgCode() == 3334  && mcGdMotherPosV0Dghter->GetPdgCode() == 3334  && mcGdMotherNegV0Dghter->GetPdgCode() == 3334)  {lAssoOmegaMinus = kTRUE; cascadeMass = 1.672; nAssoOmegaMinus++;}
             else if (mcMotherBach->GetPdgCode() == -3334 && mcGdMotherPosV0Dghter->GetPdgCode() == -3334 && mcGdMotherNegV0Dghter->GetPdgCode() == -3334) {lAssoOmegaPlus  = kTRUE; cascadeMass = 1.672; nAssoOmegaPlus++;}
             // If a proper association  exists ...	
             if(fDebug > 4){
                  cout<<"XiMinus    = "<<lAssoXiMinus   <<endl;
                  cout<<"XiPlus     = "<<lAssoXiPlus    <<endl;
                  cout<<"OmegaMinus = "<<lAssoOmegaMinus<<endl;
                  cout<<"OmegaPlus  = "<<lAssoOmegaPlus <<endl<<"----"<<endl;	
             }
             if(fDebug > 5){
                  cout<<endl;
                  cout<<"- V0 daughters - "<<endl;
                  cout<<"     + V0 Pos. / Label : "<<lblPosV0Dghter<<" - Pdg Code : "<<mcPosV0Dghter->GetTitle()<<endl;
                  cout<<"     - V0 Neg. / Label : "<<lblNegV0Dghter<<" - Pdg Code : "<<mcNegV0Dghter->GetTitle()<<endl;
                  cout<<"- Xi daughters - "<<endl;
                  cout<<"     + V0 Pos. mother / Label : "<<lblMotherPosV0Dghter<<" - Pdg Code : "<<mcMotherPosV0Dghter->GetTitle()<<endl;
                  cout<<"     - V0 Neg. mother / Label : "<<lblMotherNegV0Dghter<<" - Pdg Code : "<<mcMotherNegV0Dghter->GetTitle()<<endl;
                  cout<<"     --  Bach. / Label :"<<lblBach<<" -  Pdg Code : "<<mcBach->GetTitle()<<endl;
                  cout<<"- Xi candidate -"<<endl;
                  cout<<"    +  V0 Pos. Gd Mother / Label : "<<lblGdMotherPosV0Dghter<<" - Pdg Code : "<< mcGdMotherPosV0Dghter->GetTitle()<<endl;
                  cout<<"    -  V0 Neg. Gd Mother / Label : "<<lblGdMotherNegV0Dghter<<" - Pdg Code : "<< mcGdMotherNegV0Dghter->GetTitle()<<endl;
                  cout<<"    --  Mother Bach. / Label : "<<lblMotherBach<<" - Pdg Code    : "<<mcMotherBach->GetTitle()<<endl;
                  cout<<endl;
	     }
             lmcPt             = mcMotherBach->Pt();
             lmcRapCasc        = 0.5*TMath::Log( (mcMotherBach->Energy() + mcMotherBach->Pz()) / (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) );
             lmcEta            = mcMotherBach->Eta();
             lmcTransvRadius   = mcBach->R(); // to get the decay point of Xi, = the production vertex of Bachelor ...
             TVector3 lmcTVect3Mom( mcMotherBach->Px(), mcMotherBach->Py(), mcMotherBach->Pz() );
             lrecoPt           = xiESD->Pt();
             lrecoTransvRadius = TMath::Sqrt( xiESD->Xv() * xiESD->Xv() + xiESD->Yv() * xiESD->Yv() );
             TVector3 lrecoTVect3Mom( xiESD->Px(), xiESD->Py(), xiESD->Pz() );
             lDeltaPhiMcReco   = lmcTVect3Mom.DeltaPhi( lrecoTVect3Mom ) * 180.0/TMath::Pi();
             lmcPtPosV0Dghter  = mcPosV0Dghter->Pt() ;
             lmcPtNegV0Dghter  = mcNegV0Dghter->Pt();
             lrecoP            = xiESD->P();
             Double_t nV0mom[3] = {0. ,0. ,0. };
             Double_t pV0mom[3] = {0. ,0. ,0. };
             xiESD->GetNPxPyPz(nV0mom[0],nV0mom[1],nV0mom[2]);    
             xiESD->GetPPxPyPz(pV0mom[0],pV0mom[1],pV0mom[2]);
             lV0mom = TMath::Sqrt(TMath::Power(nV0mom[0]+pV0mom[0],2)+TMath::Power(nV0mom[1]+pV0mom[1],2)+TMath::Power(nV0mom[2]+pV0mom[2],2));
             Double_t lBachMomX = 0.; Double_t lBachMomY = 0.; Double_t lBachMomZ = 0.;
             xiESD->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
             lBachTransvMom   = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
             lnTrackTransvMom = TMath::Sqrt( nV0mom[0]*nV0mom[0] + nV0mom[1]*nV0mom[1] );
             lpTrackTransvMom = TMath::Sqrt( pV0mom[0]*pV0mom[0] + pV0mom[1]*pV0mom[1] );
            
        } else if ( fAnalysisType == "AOD" ) {

             // - Load the cascade
             xiAOD = 0x0;  xiAOD = lAODevent->GetCascade(iXi);
             if (!xiAOD) continue;
             // - Connection to daughter tracks of the current cascade
             pTrackXiaod = 0;  pTrackXiaod    = dynamic_cast<AliAODTrack*>(xiAOD->GetDaughter(0));
             nTrackXiaod = 0;  nTrackXiaod    = dynamic_cast<AliAODTrack*>(xiAOD->GetDaughter(1));
             bachTrackXiaod = 0;  bachTrackXiaod = dynamic_cast<AliAODTrack*>(xiAOD->GetDecayVertexXi()->GetDaughter(0));
             if (!pTrackXiaod || !nTrackXiaod || !bachTrackXiaod ) {  AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");  continue;  }
             lIdxPosXi = 0x0;  lIdxPosXi = (UInt_t) TMath::Abs(pTrackXiaod->GetID());
             lIdxNegXi = 0x0;  lIdxNegXi = (UInt_t) TMath::Abs(nTrackXiaod->GetID());
             lBachIdx = 0x0;  lBachIdx  = (UInt_t) TMath::Abs(bachTrackXiaod->GetID());
             // - Rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
             if(lBachIdx == lIdxNegXi) {  AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!");  continue;  }
             if(lBachIdx == lIdxPosXi) {  AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!");  continue;  }
             // - Rejection of a poor quality tracks
             if (fTrackQualityCutTPCrefit) {
                 // - Poor quality related to TPCrefit
                 if (!(pTrackXiaod->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
                 if (!(nTrackXiaod->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
                 if (!(bachTrackXiaod->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
             }
             if (fTrackQualityCutnTPCcls) {
                 lPosTPCClusters = 0;  lPosTPCClusters  = pTrackXiaod->GetTPCClusterInfo(2,1);     //->GetTPCNcls();
                 lNegTPCClusters = 0;  lNegTPCClusters  = nTrackXiaod->GetTPCClusterInfo(2,1);     //->GetTPCNcls();
                 lBachTPCClusters = 0; lBachTPCClusters = bachTrackXiaod->GetTPCClusterInfo(2,1);  //->GetTPCNcls();
                 lPosTPCFindClusters   = pTrackXiaod->GetTPCNclsF();                  // New
                 lNegTPCFindClusters   = nTrackXiaod->GetTPCNclsF();                  // New
                 lBachTPCFindClusters  = bachTrackXiaod->GetTPCNclsF();               // New
                 // - Poor quality related to TPC clusters
                 if(lPosTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Pos. track has less than 80 TPC clusters ... continue!"); continue; }
                 if(lNegTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Neg. track has less than 80 TPC clusters ... continue!"); continue; }
                 if(lBachTPCClusters < fMinnTPCcls) { AliWarning("Pb / Bach.   track has less than 80 TPC clusters ... continue!"); continue; }
                 // - Poor quality related to clusters/findable
                 if( lPosTPCFindClusters <= 0 || lNegTPCFindClusters <= 0 || lBachTPCFindClusters ) { AliWarning("Pb / Number of findable cluster <= 0 ... continue!"); continue; }
                 if ((lPosTPCClusters/lPosTPCFindClusters)    < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / V0 Pos. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
                 if ((lNegTPCClusters/lNegTPCFindClusters)    < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / V0 Neg. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
                 if ((lBachTPCClusters/lBachTPCFindClusters)  < fMinTPCcrossrawoverfindable) { AliWarning(Form("Pb / Bach. track has ratio clusters/findable < %f ... continue!",fMinTPCcrossrawoverfindable)); continue; }
             }
             etaPos  = pTrackXiaod->Eta();
             etaNeg  = nTrackXiaod->Eta();
             etaBach = bachTrackXiaod->Eta();
             // - Info over reconstructed cascades
             if( bachTrackXi->Charge() < 0 ) {
                  lInvMassXiMinus = xiAOD->MassXi();
                  lInvMassOmegaMinus = xiAOD->MassOmega();
             }
             if( bachTrackXi->Charge() >  0 ){
                  lInvMassXiPlus = xiAOD->MassXi();
                  lInvMassOmegaPlus = xiAOD->MassOmega();
             }
             lDcaXiDaughters            = xiAOD->DcaXiDaughters();
             lDcaBachToPrimVertexXi     = xiAOD->DcaBachToPrimVertex();
             lXiCosineOfPointingAngle   = xiAOD->CosPointingAngleXi( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
             lPosXi[0]                  = xiAOD->DecayVertexXiX();
             lPosXi[1]                  = xiAOD->DecayVertexXiY();
             lPosXi[2]                  = xiAOD->DecayVertexXiZ();
             lInvMassLambdaAsCascDghter = xiAOD->MassLambda();
             lDcaV0DaughtersXi          = xiAOD->DcaV0Daughters();
             lV0CosineOfPointingAngleXi = xiAOD->CosPointingAngle( lPosXi );
             lV0CosineOfPointingAngle   = xiAOD->CosPointingAngle( lBestPrimaryVtxPos );
             lPosV0Xi[0]                = xiAOD->DecayVertexV0X();
             lPosV0Xi[1]                = xiAOD->DecayVertexV0Y();
             lPosV0Xi[2]                = xiAOD->DecayVertexV0Z();
             lDcaV0ToPrimVertexXi       = xiAOD->DcaV0ToPrimVertex();
             lDcaPosToPrimVertexXi      = xiAOD->DcaPosToPrimVertex();
             lDcaNegToPrimVertexXi      = xiAOD->DcaNegToPrimVertex();
             lChargeXi                  = xiAOD->ChargeXi();
             // - PID Information
                // Combined VO-positive-daughter PID
                // Combined bachelor PID
                /* 
                AliPID bachPidXi;       bachPidXi.SetPriors(    lPriorsGuessXi    );
                AliPID bachPidOmega;    bachPidOmega.SetPriors( lPriorsGuessOmega );
                if ( bachTrackXi->IsOn(AliESDtrack::kESDpid) ) {  // Combined PID exists
                     Double_t r[10] = {0.}; bachTrackXi->GetESDpid(r);
                     bachPidXi.SetProbabilities(r);
                     bachPidOmega.SetProbabilities(r);
                     // Check if the bachelor track is a pion
                     ppionBach = bachPidXi.GetProbability(AliPID::kPion);
                     if (ppionBach > bachPidXi.GetProbability(AliPID::kElectron) &&
                         ppionBach > bachPidXi.GetProbability(AliPID::kMuon)     &&
                         ppionBach > bachPidXi.GetProbability(AliPID::kKaon)     &&
                         ppionBach > bachPidXi.GetProbability(AliPID::kProton)   )     lIsBachelorPion = kTRUE;
                     // Check if the bachelor track is a kaon
                     pkaonBach = bachPidOmega.GetProbability(AliPID::kKaon);
                     if (pkaonBach > bachPidOmega.GetProbability(AliPID::kElectron) &&
                         pkaonBach > bachPidOmega.GetProbability(AliPID::kMuon)     &&
                         pkaonBach > bachPidOmega.GetProbability(AliPID::kPion)     &&
                         pkaonBach > bachPidOmega.GetProbability(AliPID::kProton)   )  lIsBachelorKaon = kTRUE;
                }// end if bachelor track with existing combined PID
                */
                // - TPC PID: 4-sigma bands on Bethe-Bloch curve
                // Bachelor
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXiaod,AliPID::kKaon)) < fTPCPIDsigma) lIsBachelorKaonForTPC = kTRUE;
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXiaod,AliPID::kPion)) < fTPCPIDsigma) lIsBachelorPionForTPC = kTRUE;
                // Negative V0 daughter
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXiaod,AliPID::kPion   )) < fTPCPIDsigma) lIsNegPionForTPC   = kTRUE;
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXiaod,AliPID::kProton )) < fTPCPIDsigma) lIsNegProtonForTPC = kTRUE;
                // Positive V0 daughter
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXiaod,AliPID::kPion   )) < fTPCPIDsigma) lIsPosPionForTPC   = kTRUE;
                if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXiaod,AliPID::kProton )) < fTPCPIDsigma) lIsPosProtonForTPC = kTRUE;
                /*
                const AliExternalTrackParam *pInnerWallTrackXi    = pTrackXiaod    ->GetInnerParam(); // Do not use GetTPCInnerWall
                const AliExternalTrackParam *nInnerWallTrackXi    = nTrackXiaod    ->GetInnerParam();
                const AliExternalTrackParam *bachInnerWallTrackXi = bachTrackXiaod ->GetInnerParam();
                if(pInnerWallTrackXi && nInnerWallTrackXi && bachInnerWallTrackXi ){
                    Double_t pMomInnerWall    = pInnerWallTrackXi   ->GetP();
                    Double_t nMomInnerWall    = nInnerWallTrackXi   ->GetP();
                    Double_t bachMomInnerWall = bachInnerWallTrackXi->GetP();
                    // Bachelor
                    if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXiaod,AliPID::kPion)) < 3)                              lIsBachelorPionForTPC = kTRUE;
                    if (bachMomInnerWall < 0.350  && TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXiaod,AliPID::kKaon)) < 5) lIsBachelorKaonForTPC = kTRUE;
                    if (bachMomInnerWall > 0.350  && TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXiaod,AliPID::kKaon)) < 3) lIsBachelorKaonForTPC = kTRUE;
                    // Negative V0 daughter
                    if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXiaod,AliPID::kPion   )) < 3  )                           lIsNegPionForTPC   = kTRUE;
                    if (nMomInnerWall < 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXiaod,AliPID::kProton ) ) < 5 )   lIsNegProtonForTPC = kTRUE;
                    if (nMomInnerWall > 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXiaod,AliPID::kProton ) ) < 3 )   lIsNegProtonForTPC = kTRUE;
                    // Positive V0 daughter
                    if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXiaod,AliPID::kPion   )) < 3 )                            lIsPosPionForTPC   = kTRUE;
                    if (pMomInnerWall < 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXiaod,AliPID::kProton )) < 5)     lIsPosProtonForTPC = kTRUE;
                    if (pMomInnerWall > 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXiaod,AliPID::kProton )) < 3)     lIsPosProtonForTPC = kTRUE;
                }*/
                // - PID proba Vs Pt(Bach)
                lblBachForPID = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
                mcBachForPIDaod = 0x0; mcBachForPIDaod   = (AliAODMCParticle*) arrayMC->At( lblBachForPID );
                lmcPtBach = mcBachForPIDaod->Pt();
                // - MC perfect PID
                if (mcBachForPIDaod->PdgCode() == -211) lIsBachelorMCPiMinus = kTRUE;
                if (mcBachForPIDaod->PdgCode() ==  211) lIsBachelorMCPiPlus  = kTRUE;
                if (mcBachForPIDaod->PdgCode() == -321) lIsBachelorMCKMinus  = kTRUE;
                if (mcBachForPIDaod->PdgCode() ==  321) lIsBachelorMCKPlus   = kTRUE;
             // - MC association (care : lots of "continue;" below this line)
             if (fDebug > 5) cout<<"MC EventNumber : "<<lMCevent->Header()->GetEvent()<<" / MC event Number in Run : "<<lMCevent->Header()->GetEventNrInRun()<<endl;
             // - Level of the V0 daughters
             lblPosV0Dghter = 0; lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXiaod->GetLabel() );
             lblNegV0Dghter = 0; lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXiaod->GetLabel() );
             mcPosV0Dghteraod = 0x0; mcPosV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblPosV0Dghter );
             mcNegV0Dghteraod = 0x0; mcNegV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblNegV0Dghter );
             // - Level of the Xi daughters
             lblMotherPosV0Dghter = 0;  lblMotherPosV0Dghter = mcPosV0Dghteraod->GetMother();   
             lblMotherNegV0Dghter = 0;  lblMotherNegV0Dghter = mcNegV0Dghteraod->GetMother();
             if (lblMotherPosV0Dghter != lblMotherNegV0Dghter) continue; // same mother
             if (lblMotherPosV0Dghter < 0 ) continue;                    // this particle is primary, no mother
             if (lblMotherNegV0Dghter < 0 ) continue;                    // this particle is primary, no mother
                                                                         // mothers = Lambda candidate ... a priori
             mcMotherPosV0Dghteraod = 0x0; mcMotherPosV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblMotherPosV0Dghter );
             mcMotherNegV0Dghteraod = 0x0; mcMotherNegV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblMotherNegV0Dghter );  
             lblBach = 0;  lblBach = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
             mcBachaod = 0x0; mcBachaod = (AliAODMCParticle*) arrayMC->At( lblBach );
             // - Level of Xi candidate
             lblGdMotherPosV0Dghter = 0; lblGdMotherPosV0Dghter = mcMotherPosV0Dghteraod->GetMother() ;
             lblGdMotherNegV0Dghter = 0; lblGdMotherNegV0Dghter = mcMotherNegV0Dghteraod->GetMother() ;
             if (lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter ) continue;
             if (lblGdMotherPosV0Dghter < 0 ) continue;                    // primary lambda ...
             if (lblGdMotherNegV0Dghter < 0 ) continue;                    // primary lambda ...
                                                                           // Gd mothers = Xi candidate ... a priori
             mcGdMotherPosV0Dghteraod = 0x0; mcGdMotherPosV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblGdMotherPosV0Dghter );
             mcGdMotherNegV0Dghteraod = 0x0; mcGdMotherNegV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblGdMotherNegV0Dghter );
             lblMotherBach = 0; lblMotherBach = (Int_t) TMath::Abs( mcBachaod->GetMother() );
             if (lblMotherBach != lblGdMotherPosV0Dghter ) continue; //same mother for bach and V0 daughters
             mcMotherBachaod = 0x0; mcMotherBachaod = (AliAODMCParticle*) arrayMC->At( lblMotherBach );
             // - Check if cascade is primary
             if (!(mcMotherBachaod->IsPhysicalPrimary())) continue;
             // - Manage boolean for association
             if      (mcMotherBachaod->GetPdgCode() == 3312  && mcGdMotherPosV0Dghteraod->GetPdgCode() == 3312  && mcGdMotherNegV0Dghteraod->GetPdgCode() == 3312 ) {lAssoXiMinus = kTRUE;    cascadeMass = 1.321; nAssoXiMinus++; }
             else if (mcMotherBachaod->GetPdgCode() == -3312 && mcGdMotherPosV0Dghteraod->GetPdgCode() == -3312 && mcGdMotherNegV0Dghteraod->GetPdgCode() == -3312) {lAssoXiPlus = kTRUE;     cascadeMass = 1.321; nAssoXiPlus++; }
             else if (mcMotherBachaod->GetPdgCode() == 3334  && mcGdMotherPosV0Dghteraod->GetPdgCode() == 3334  && mcGdMotherNegV0Dghteraod->GetPdgCode() == 3334 ) {lAssoOmegaMinus = kTRUE; cascadeMass = 1.672; nAssoOmegaMinus++; }
             else if (mcMotherBachaod->GetPdgCode() == -3334 && mcGdMotherPosV0Dghteraod->GetPdgCode() == -3334 && mcGdMotherNegV0Dghteraod->GetPdgCode() == -3334) {lAssoOmegaPlus = kTRUE;  cascadeMass = 1.672; nAssoOmegaPlus++; }
             lmcPt              = mcMotherBach->Pt();
             lmcRapCasc         = 0.5*TMath::Log( (mcMotherBach->Energy() + mcMotherBach->Pz()) / (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) );
             lmcEta             = mcMotherBach->Eta();
             decayCascX = 0.;  decayCascX = mcBachaod->Xv();
             decayCascY = 0.;  decayCascY = mcBachaod->Yv();
             lmcTransvRadius    = TMath::Sqrt(decayCascX*decayCascX+decayCascY*decayCascY); // decay point of Xi, = the production vertex of Bachelor ...
             TVector3 lmcTVect3Mom( mcMotherBach->Px(), mcMotherBach->Py(), mcMotherBach->Pz() );
             Double_t xiMomX    = xiAOD->MomXiX();
             Double_t xiMomY    = xiAOD->MomXiY();
             Double_t xiMomZ    = xiAOD->MomXiZ();
             lrecoPt            = TMath::Sqrt( xiMomX*xiMomX   + xiMomY*xiMomY ); 
             lrecoTransvRadius  = TMath::Sqrt( xiAOD->DecayVertexXiX() * xiAOD->DecayVertexXiX() + xiAOD->DecayVertexXiY() * xiAOD->DecayVertexXiY() );
             TVector3 lrecoTVect3Mom( xiMomX, xiMomY, xiMomZ );
             lDeltaPhiMcReco    = lmcTVect3Mom.DeltaPhi( lrecoTVect3Mom ) * 180.0/TMath::Pi();
             lmcPtPosV0Dghter   = mcPosV0Dghteraod->Pt() ;
             lmcPtNegV0Dghter   = mcNegV0Dghteraod->Pt();
             lrecoP             = TMath::Sqrt( xiMomX*xiMomX   + xiMomY*xiMomY   + xiMomZ*xiMomZ );;
             lV0momX = 0.; lV0momX   = xiAOD->MomV0X();
             lV0momY = 0.; lV0momY   = xiAOD->MomV0Y();
             lV0momZ = 0.; lV0momZ   = xiAOD->MomV0Z();
             lV0mom = TMath::Sqrt(TMath::Power(lV0momX,2)+TMath::Power(lV0momY,2)+TMath::Power(lV0momZ,2));
             lBachMomX = 0.; lBachMomX = xiAOD->MomBachX();
             lBachMomY = 0.; lBachMomY = xiAOD->MomBachY();
             lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
             lV0NMomX = 0.; lV0NMomX = xiAOD->MomNegX();
             lV0NMomY = 0.; lV0NMomY = xiAOD->MomNegY();
             lV0PMomX = 0.; lV0PMomX = xiAOD->MomPosX();
             lV0PMomY = 0.; lV0PMomY = xiAOD->MomPosY();
             lnTrackTransvMom = TMath::Sqrt( lV0NMomX*lV0NMomX   + lV0NMomY*lV0NMomY );
             lpTrackTransvMom = TMath::Sqrt( lV0PMomX*lV0PMomX   + lV0PMomY*lV0PMomY );
            
        }
        lXiRadius   = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
        lV0RadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] ); 
        // - Cut on pt of the three daughter tracks
        if (lBachTransvMom   < fMinPtCutOnDaughterTracks) continue;
        if (lpTrackTransvMom < fMinPtCutOnDaughterTracks) continue;
        if (lnTrackTransvMom < fMinPtCutOnDaughterTracks) continue;
        // - Cut on pseudorapidity of the three daughter tracks
        if (TMath::Abs(etaBach) > fEtaCutOnDaughterTracks) continue;
        if (TMath::Abs(etaPos)  > fEtaCutOnDaughterTracks) continue;
        if (TMath::Abs(etaNeg)  > fEtaCutOnDaughterTracks) continue;
        // - Extra-selection for cascade candidates
        if (fkExtraSelections) {
                if (lDcaXiDaughters > 0.3) continue;                                              // in AliCascadeVertexer
                if (lXiCosineOfPointingAngle < 0.999 ) continue;                                  // in AliCascadeVertexer
                if (lDcaV0ToPrimVertexXi < 0.05) continue;                                        // in AliCascadeVertexer
                if (lDcaBachToPrimVertexXi < 0.03) continue;                                      // in AliCascadeVertexer
                if (lDcaV0DaughtersXi > 1.) continue;                                             // in AliV0vertexer
                if ((fCollidingSystem == 0) && (lV0CosineOfPointingAngleXi < 0.998)) continue; // in AliV0vertexer
                if ((fCollidingSystem == 1) && (lV0CosineOfPointingAngle < 0.998)) continue;  // in AliV0vertexer
                if (lDcaPosToPrimVertexXi < 0.1) continue;                                        // in AliV0vertexer
                if (lDcaNegToPrimVertexXi < 0.1) continue;                                        // in AliV0vertexer
                if (lXiRadius < .9) continue;                                                     // in AliCascadeVertexer
                if (lV0RadiusXi < 0.9) continue;                                                  // in AliV0vertexer
        }
        // - Fill combined PID TH1s
        if (lIsBachelorPion)   f2dHistPIDprobaPionVsMCPtBach->Fill( lmcPtBach, ppionBach );
        if (lIsBachelorKaon)   f2dHistPIDprobaKaonVsMCPtBach->Fill( lmcPtBach, pkaonBach );
        // - No association, skip the rest of the code
        if(!lAssoXiMinus && !lAssoXiPlus && !lAssoOmegaMinus && !lAssoOmegaPlus) continue; 
        // - Proper time         
        // -- For cascade (reconstructed)   
        lctau = 0.;
        lctau = TMath::Sqrt(TMath::Power((lPosXi[0]-lBestPrimaryVtxPos[0]),2)+TMath::Power((lPosXi[1]-lBestPrimaryVtxPos[1]),2)+TMath::Power((lPosXi[2]-lBestPrimaryVtxPos[2]),2));
        if (lrecoP!=0) lctau = lctau*cascadeMass/lrecoP;   
        else           lctau = -1.;
        // -- For Lambda (reconstructed)
        distV0Xi = 0.; distV0Xi = TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2)+TMath::Power((lPosV0Xi[2]-lPosXi[2]),2)); 
        lctauV0 = -1.; if (lV0mom!=0) lctauV0 = distV0Xi*lambdaMass/lV0mom; 
        // - Distance of the V0 and cascade vertices
        distTV0Xi = 0.; distTV0Xi = TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2));
        // - Fill histos for the cascade candidates associated with MC
        if (lChargeXi < 0 && lAssoXiMinus){	
                fHistAsMCMassXiMinus->Fill(lInvMassXiMinus);
                if (lIsBachelorPion) f2dHistAsMCandCombPIDGenPtVsGenYXiMinus->Fill(lmcPt, lmcRapCasc);
                f2dHistAsMCGenPtVsGenYXiMinus->Fill(lmcPt, lmcRapCasc);
                fHistAsMCGenEtaXiMinus->Fill(lmcEta);
                f2dHistAsMCResPtXiMinus->Fill(lmcPt, (lrecoPt - lmcPt) / lmcPt);
                f2dHistAsMCResRXiMinus->Fill(lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius) / lmcTransvRadius);
                f2dHistAsMCResPhiXiMinus->Fill(lmcPt, lDeltaPhiMcReco);
                f2dHistAsMCptProtonMCptXiMinus->Fill(lmcPt, lmcPtPosV0Dghter);
        }	
        else if (lChargeXi > 0 && lAssoXiPlus){	
                fHistAsMCMassXiPlus->Fill(lInvMassXiPlus);
                if (lIsBachelorPion) f2dHistAsMCandCombPIDGenPtVsGenYXiPlus->Fill(lmcPt, lmcRapCasc);
                f2dHistAsMCGenPtVsGenYXiPlus->Fill(lmcPt, lmcRapCasc);
                fHistAsMCGenEtaXiPlus->Fill(lmcEta);
                f2dHistAsMCResPtXiPlus->Fill(lmcPt, (lrecoPt - lmcPt) / lmcPt);
                f2dHistAsMCResRXiPlus->Fill(lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius) / lmcTransvRadius);
                f2dHistAsMCResPhiXiPlus->Fill(lmcPt, lDeltaPhiMcReco);
                f2dHistAsMCptAntiprotonMCptXiPlus->Fill(lmcPt, lmcPtNegV0Dghter);
        }
        else if( lChargeXi < 0 && lAssoOmegaMinus){	
                fHistAsMCMassOmegaMinus->Fill(lInvMassOmegaMinus);
                if (lIsBachelorKaon) f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus->Fill(lmcPt, lmcRapCasc);
                f2dHistAsMCGenPtVsGenYOmegaMinus->Fill(lmcPt, lmcRapCasc);
                fHistAsMCGenEtaOmegaMinus->Fill(lmcEta);
                f2dHistAsMCResPtOmegaMinus->Fill(lmcPt, (lrecoPt - lmcPt) / lmcPt );
                f2dHistAsMCResROmegaMinus->Fill(lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius) / lmcTransvRadius);
                f2dHistAsMCResPhiOmegaMinus->Fill( lmcPt, lDeltaPhiMcReco);
                f2dHistAsMCptProtonMCptOmegaMinus->Fill(lmcPt, lmcPtPosV0Dghter);
        }	
        else if( lChargeXi > 0 && lAssoOmegaPlus){	
                fHistAsMCMassOmegaPlus->Fill(lInvMassOmegaPlus);
                if (lIsBachelorKaon) f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus->Fill(lmcPt, lmcRapCasc);
                f2dHistAsMCGenPtVsGenYOmegaPlus->Fill(lmcPt, lmcRapCasc);
                fHistAsMCGenEtaOmegaPlus->Fill(lmcEta);
                f2dHistAsMCResPtOmegaPlus->Fill(lmcPt, (lrecoPt - lmcPt) / lmcPt);
                f2dHistAsMCResROmegaPlus->Fill(lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius) / lmcTransvRadius);
                f2dHistAsMCResPhiOmegaPlus->Fill(lmcPt, lDeltaPhiMcReco);
                f2dHistAsMCptAntiprotonMCptOmegaPlus->Fill(lmcPt, lmcPtNegV0Dghter);
        }
        fHistPosV0TPCClusters->Fill( lPosTPCClusters );
        fHistNegV0TPCClusters->Fill( lNegTPCClusters );
        fHistBachTPCClusters->Fill( lBachTPCClusters );
        // - Fill containers
        // - Filling the AliCFContainer (optimisation of topological selections + systematics)
        Double_t lContainerCutVars[20] = {0.0};
        lContainerCutVars[0]  = lDcaXiDaughters;
        lContainerCutVars[1]  = lDcaBachToPrimVertexXi;
        lContainerCutVars[2]  = lXiCosineOfPointingAngle;
        lContainerCutVars[3]  = lXiRadius;
        lContainerCutVars[4]  = lInvMassLambdaAsCascDghter;
        lContainerCutVars[5]  = lDcaV0DaughtersXi;
        if      (fCollidingSystem == 0) lContainerCutVars[6]  = lV0CosineOfPointingAngleXi;
        else if (fCollidingSystem == 1) lContainerCutVars[6]  = lV0CosineOfPointingAngle;
        lContainerCutVars[7]  = lV0RadiusXi;
        lContainerCutVars[8]  = lDcaV0ToPrimVertexXi;	
        lContainerCutVars[9]  = lDcaPosToPrimVertexXi;
        lContainerCutVars[10] = lDcaNegToPrimVertexXi;
        lContainerCutVars[13] = lmcPt;
        lContainerCutVars[16] = lctau;
        lContainerCutVars[17] = lctauV0;
        lContainerCutVars[18] = distTV0Xi;
        // All cases should be covered below
        if (lChargeXi < 0 && lAssoXiMinus) {
                     lContainerCutVars[11] = lInvMassXiMinus;
                     lContainerCutVars[12] = lInvMassOmegaMinus;
                     lContainerCutVars[14] = lmcRapCasc;
                     lContainerCutVars[15] = -1.;
                     lContainerCutVars[19] = TMath::Abs(lInnerWallMomCascDghters[0]);
                     if (lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) fCFContAsCascadeCuts->Fill(lContainerCutVars,0);
        }
        if (lChargeXi > 0 && lAssoXiPlus) {
                     lContainerCutVars[11] = lInvMassXiPlus;
                     lContainerCutVars[12] = lInvMassOmegaPlus;
                     lContainerCutVars[14] = lmcRapCasc;
                     lContainerCutVars[15] = -1.; 
                     lContainerCutVars[19] = TMath::Abs(lInnerWallMomCascDghters[1]);
                     if (lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC) fCFContAsCascadeCuts->Fill(lContainerCutVars,1);
        }
        if (lChargeXi < 0 && lAssoOmegaMinus) {
                     lContainerCutVars[11] = lInvMassXiMinus;
                     lContainerCutVars[12] = lInvMassOmegaMinus;
                     lContainerCutVars[14] = -1.;
                     lContainerCutVars[15] = lmcRapCasc;
                     lContainerCutVars[19] = TMath::Abs(lInnerWallMomCascDghters[0]);
                     if (lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) fCFContAsCascadeCuts->Fill(lContainerCutVars,2);
        }
        if(lChargeXi > 0 && lAssoOmegaPlus) {
                     lContainerCutVars[11] = lInvMassXiPlus;
                     lContainerCutVars[12] = lInvMassOmegaPlus;
                     lContainerCutVars[14] = -1.;
                     lContainerCutVars[15] = lmcRapCasc;
                     lContainerCutVars[19] = TMath::Abs(lInnerWallMomCascDghters[1]);
                     if (lIsBachelorKaonForTPC && lIsNegProtonForTPC && lIsPosPionForTPC) fCFContAsCascadeCuts->Fill(lContainerCutVars,3);
        }
        // - Filling the AliCFContainers related to PID
        Double_t lContainerPIDVars[3] = {0.0};
        // Xi Minus		
        if (lChargeXi < 0 && lAssoXiMinus) {
              lContainerPIDVars[0] = lmcPt;
              lContainerPIDVars[1] = lInvMassXiMinus;
              lContainerPIDVars[2] = lmcRapCasc;
              fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 0);                                                                      // No PID
              if (lIsBachelorPionForTPC) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 1);                                           // TPC PID / 4-#sigma cut on Bachelor track
              if (lIsBachelorPionForTPC && lIsPosProtonForTPC) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 2);                     // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks	
              if (lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
              if (lIsBachelorPion) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 4);                                                 // Comb. PID / Bachelor
              if (lIsBachelorPion && lIsPosInXiProton) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 5);                             // Comb. PID / Bachelor+Baryon
              if (lIsBachelorPion && lIsPosInXiProton && lIsNegInXiPion) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 6);           // Comb. PID / Bachelor+Baryon+Meson
        }	
        lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; 
        // Xi Plus		
        if( lChargeXi > 0 && lAssoXiPlus ) {
              lContainerPIDVars[0] = lmcPt;
              lContainerPIDVars[1] = lInvMassXiPlus;
    	      lContainerPIDVars[2] = lmcRapCasc;
              fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 0);                                                                      // No PID
              if (lIsBachelorPionForTPC) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 1);                                           // TPC PID / 4-#sigma cut on Bachelor track
              if (lIsBachelorPionForTPC && lIsNegProtonForTPC) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 2);                     // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
              if (lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
              if (lIsBachelorPion ) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 4);                                                // Comb. PID / Bachelor
              if (lIsBachelorPion && lIsNegInXiProton) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 5);                             // Comb. PID / Bachelor+Baryon
              if (lIsBachelorPion && lIsNegInXiProton && lIsPosInXiPion) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 6);           // Comb. PID / Bachelor+Baryon+Meson
        }	
        lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; 
        // Omega Minus		
        if( lChargeXi < 0 && lAssoOmegaMinus ) {
              lContainerPIDVars[0] = lmcPt;
              lContainerPIDVars[1] = lInvMassOmegaMinus;
              lContainerPIDVars[2] = lmcRapCasc;		
              fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 0);                                                                      // No PID
              if (lIsBachelorKaonForTPC) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 1);                                           // TPC PID / 4-#sigma cut on Bachelor track
              if (lIsBachelorKaonForTPC && lIsPosProtonForTPC) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 2);                     // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
              if (lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
              if (lIsBachelorKaon) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 4);                                                 // Comb. PID / Bachelor
              if (lIsBachelorKaon && lIsPosInOmegaProton) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 5);                          // Comb. PID / Bachelor+Baryon
              if (lIsBachelorKaon && lIsPosInOmegaProton && lIsNegInOmegaPion) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 6);     // Comb. PID / Bachelor+Baryon+Meson
        }	
        lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; 
        // Omega Plus		
        if( lChargeXi > 0 && lAssoOmegaPlus) {
              lContainerPIDVars[0] = lmcPt;
              lContainerPIDVars[1] = lInvMassOmegaPlus;
              lContainerPIDVars[2] = lmcRapCasc;		
              fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 0);                                                                      // No PID
              if (lIsBachelorKaonForTPC) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 1);                                           // TPC PID / 4-#sigma cut on Bachelor track
              if (lIsBachelorKaonForTPC && lIsNegProtonForTPC) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 2);                     // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
              if (lIsBachelorKaonForTPC && lIsNegProtonForTPC && lIsPosPionForTPC) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
              if (lIsBachelorKaon) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 4);                                                 // Comb. PID / Bachelor
              if (lIsBachelorKaon && lIsNegInOmegaProton) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 5);                          // Comb. PID / Bachelor+Baryon
              if (lIsBachelorKaon && lIsNegInOmegaProton && lIsPosInOmegaPion) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 6);     // Comb. PID / Bachelor+Baryon+Meson
        }	
	
   }// End of loop over reconstructed cascades
 
   fHistnAssoXiMinus->Fill(nAssoXiMinus);
   fHistnAssoXiPlus->Fill(nAssoXiPlus);
   fHistnAssoOmegaMinus->Fill(nAssoOmegaMinus);
   fHistnAssoOmegaPlus->Fill(nAssoOmegaPlus);  
   // Post output data.
   PostData(1, fListHistCascade);
   PostData(2, fCFContCascadePIDAsXiMinus);
   PostData(3, fCFContCascadePIDAsXiPlus);
   PostData(4, fCFContCascadePIDAsOmegaMinus);
   PostData(5, fCFContCascadePIDAsOmegaPlus);
   PostData(6, fCFContAsCascadeCuts);

}      


//________________________________________________________________________
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// - TERMINATE
void AliAnalysisTaskCheckPerformanceCascadepp::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
	
 /* TList *cRetrievedList = 0x0;
  cRetrievedList = (TList*)GetOutputData(1);
  if(!cRetrievedList) {
	Printf("ERROR - AliAnalysisTaskCheckPerformanceCascadepp : ouput data container list not available\n");
	return;
  }	
	
   
  TCanvas *canCheckPerformanceCascade = new TCanvas("AliAnalysisTaskCheckPerformanceCascadepp","Multiplicity",10,10,510,510);
  canCheckPerformanceCascade->cd(1)->SetLogy();

 */
}
