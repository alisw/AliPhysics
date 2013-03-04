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
//
//
//              Adapted to pp 2.76 analysis: D. Colella, domenico.colella@ba.infn.it (Nov. 2012):
//                        - added new and removed other histograms 
//                        - Physics selection moved here (mainly for normalization in the efficiency calcuation)
//                        - Centrality selection deleted
//                        - 3DHisto denominator moved before any event selection for Normalization
//                        - injected and natural part of MC selection removed
// 
//
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
//#include "AliV0vertexer.h"
//#include "AliCascadeVertexer.h"
#include "AliESDEvent.h"
#include "AliESDcascade.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h" 
#include "AliAnalysisTaskCheckPerformanceCascadepp276.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskCheckPerformanceCascadepp276)



//________________________________________________________________________________________
AliAnalysisTaskCheckPerformanceCascadepp276::AliAnalysisTaskCheckPerformanceCascadepp276() 
: AliAnalysisTaskSE(), // <- take care to AliAnalysisTask( empty )
    fAnalysisType                  ("ESD"), 
    fESDtrackCuts                  (0), 
    fPIDResponse                   (0),
    fkRerunV0CascVertexers         (0),
    fkSDDselectionOn               (kTRUE),
    fkQualityCutZprimVtxPos        (kTRUE),
    fkRejectEventPileUp            (kTRUE),
    fkQualityCutNoTPConlyPrimVtx   (kTRUE),
    fkQualityCutTPCrefit           (kTRUE),
    fkQualityCutnTPCcls            (kTRUE),
    fwithSDD                       (kTRUE),
    fMinnTPCcls                    (0),
    fkExtraSelections              (0),
    fVtxRange                      (0),
    fVtxRangeMin                   (0),
    fApplyAccCut                   (0),
    fMinPtCutOnDaughterTracks      (0),
    fEtaCutOnDaughterTracks        (0),
    
    // - Plots initialisation
    fListHistCascade(0),

    // - General Plots
    // Cascade multiplicity plots
    fHistCascadeMultiplicityBeforeAnySel(0),
    fHistCascadeMultiplicityAfterSDDSel(0),
    fHistCascadeMultiplicityAfterPhysicsSel(0),
    fHistCascadeMultiplicityForSelEvtNoTPCOnly(0),
    fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup(0),
    fHistCascadeMultiplicityAfterVertexCutSel(0),
    fHistnXiPlusPerEvTot(0),                  // After any event selections, in all the eta and pt range
    fHistnXiMinusPerEvTot(0),                 // After any event selections, in all the eta and pt range
    fHistnOmegaPlusPerEvTot(0),               // After any event selections, in all the eta and pt range
    fHistnOmegaMinusPerEvTot(0),              // After any event selections, in all the eta and pt range
    fHistnXiPlusPerEv(0),                     // After any event selections, in the detector acceptance and over a pt minimum
    fHistnXiMinusPerEv(0),                    // After any event selections, in the detector acceptance and over a pt minimum
    fHistnOmegaPlusPerEv(0),                  // After any event selections, in the detector acceptance and over a pt minimum
    fHistnOmegaMinusPerEv(0),                 // After any event selections, in the detector acceptance and over a pt minimum
    fHistnAssoXiMinus(0),                     // For the Reconstructed-Associated cascades 
    fHistnAssoXiPlus(0),                      // For the Reconstructed-Associated cascades 
    fHistnAssoOmegaMinus(0),                  // For the Reconstructed-Associated cascades 
    fHistnAssoOmegaPlus(0),                   // For the Reconstructed-Associated cascades 
    // Tracks multiplicity plots
    fHistTrackMultiplicityBeforeAnySel(0),
    fHistTrackMultiplicityAfterSDDSel(0),
    fHistTrackMultiplicityAfterPhysicsSel(0),
    fHistTrackMultiplicityForSelEvtNoTPCOnly(0),
    fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup(0),
    fHistTrackMultiplicityAfterVertexCutSel(0),
    // Vertex position plots (BestVertex)
    fHistPVx(0),                              // After any selections but before |Z| < 10 cm
    fHistPVy(0),                              // After any selections but before |Z| < 10 cm
    fHistPVz(0),                              // After any selections but before |Z| < 10 cm
    fHistPVxAnalysis(0),                      // After any event selections
    fHistPVyAnalysis(0),                      // After any event selections
    fHistPVzAnalysis(0),                      // After any event selections
    // - Plots before Physics Selection
    f3dHistGenPtVsGenYvsNtracksXiMinus(0),    // After the SDD event selection (For efficinecy calculation)
    f3dHistGenPtVsGenctauvsYXiMinus(0),       // After the SDD event selection (For efficinecy calculation)
    f3dHistGenPtVsGenYvsNtracksXiPlus(0),     // After the SDD event selection (For efficinecy calculation)
    f3dHistGenPtVsGenctauvsYXiPlus(0),        // After the SDD event selection (For efficinecy calculation)
    f3dHistGenPtVsGenYvsNtracksOmegaMinus(0), // After the SDD event selection (For efficinecy calculation)
    f3dHistGenPtVsGenctauvsYOmegaMinus(0),    // After the SDD event selection (For efficinecy calculation)
    f3dHistGenPtVsGenYvsNtracksOmegaPlus(0),  // After the SDD event selection (For efficinecy calculation)
    f3dHistGenPtVsGenctauvsYOmegaPlus(0),     // After the SDD event selection (For efficinecy calculation)
    // - Generated cascade plots
    // After all the event selections 
    //Xi-
    fHistEtaGenCascXiMinus(0),                // In all the eta and pt range (as they are generated)
    fHistThetaGenCascXiMinus(0),              // In all the eta and pt range (as they are generated)
    f3dHistGenPtVsGenYvsNtracksXiMinusPhysEff(0),    // 
    f3dHistGenPtVsGenctauvsYXiMinusPhysEff(0),       // 
    f2dHistGenPtVsGenYFdblXiMinus(0),         // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaLambdaXiMinus(0),               // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBachXiMinus(0),                 // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaMesDghterXiMinus(0),            // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBarDghterXiMinus(0),            // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBachXiMinus(0),                    // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtMesDghterXiMinus(0),               // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBarDghterXiMinus(0),               // In the detector acceptance and over a pt minimum (Findable particle)
    //Xi+
    fHistEtaGenCascXiPlus(0),                 // In all the eta and pt range (as they are generated)
    fHistThetaGenCascXiPlus(0),               // In all the eta and pt range (as they are generated)
    f3dHistGenPtVsGenYvsNtracksXiPlusPhysEff(0),    // 
    f3dHistGenPtVsGenctauvsYXiPlusPhysEff(0),       //
    f2dHistGenPtVsGenYFdblXiPlus(0),          // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaLambdaXiPlus(0),                // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBachXiPlus(0),                  // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaMesDghterXiPlus(0),             // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBarDghterXiPlus(0),             // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBachXiPlus(0),                     // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtMesDghterXiPlus(0),                // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBarDghterXiPlus(0),                // In the detector acceptance and over a pt minimum (Findable particle)
    //Omega-
    fHistEtaGenCascOmegaMinus(0),             // In all the eta and pt range (as they are generated)
    fHistThetaGenCascOmegaMinus(0),           // In all the eta and pt range (as they are generated)
    f3dHistGenPtVsGenYvsNtracksOmegaMinusPhysEff(0),    // 
    f3dHistGenPtVsGenctauvsYOmegaMinusPhysEff(0),       //
    f2dHistGenPtVsGenYFdblOmegaMinus(0),      // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaLambdaOmegaMinus(0),            // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBachOmegaMinus(0),              // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaMesDghterOmegaMinus(0),         // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBarDghterOmegaMinus(0),         // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBachOmegaMinus(0),                 // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtMesDghterOmegaMinus(0),            // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBarDghterOmegaMinus(0),            // In the detector acceptance and over a pt minimum (Findable particle)
    //Omega+      
    fHistEtaGenCascOmegaPlus(0),              // In all the eta and pt range (as they are generated)
    fHistThetaGenCascOmegaPlus(0),            // In all the eta and pt range (as they are generated)
    f3dHistGenPtVsGenYvsNtracksOmegaPlusPhysEff(0),    // 
    f3dHistGenPtVsGenctauvsYOmegaPlusPhysEff(0),       //
    f2dHistGenPtVsGenYFdblOmegaPlus(0),       // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaLambdaOmegaPlus(0),             // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBachOmegaPlus(0),               // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaMesDghterOmegaPlus(0),          // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBarDghterOmegaPlus(0),          // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBachOmegaPlus(0),                  // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtMesDghterOmegaPlus(0),             // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBarDghterOmegaPlus(0),             // In the detector acceptance and over a pt minimum (Findable particle)

    // - Associated to MC cascade plots
    fHistMassXiMinus(0),                      // For the Reconstructed-Associated cascades
    fHistMassXiPlus(0),                       // For the Reconstructed-Associated cascades
    fHistMassOmegaMinus(0),                   // For the Reconstructed-Associated cascades
    fHistMassOmegaPlus(0),                    // For the Reconstructed-Associated cascades
    // Effective mass histos with combined PID
    fHistMassWithCombPIDXiMinus(0),           
    fHistMassWithCombPIDXiPlus(0),
    fHistMassWithCombPIDOmegaMinus(0), 
    fHistMassWithCombPIDOmegaPlus(0),	
    // PID Probability versus MC Pt(bachelor track)
    f2dHistPIDprobaKaonVsMCPtBach(0), f2dHistPIDprobaPionVsMCPtBach(0),
    // Effective mass histos with perfect MC PID on the bachelor
    fHistMassWithMcPIDXiMinus(0), fHistMassWithMcPIDXiPlus(0),
    fHistMassWithMcPIDOmegaMinus(0), fHistMassWithMcPIDOmegaPlus(0),
    // Effective mass histos for the cascade candidates associated with MC
    fHistAsMCMassXiMinus(0),		
    fHistAsMCMassXiPlus(0),		
    fHistAsMCMassOmegaMinus(0),
    fHistAsMCMassOmegaPlus(0),
    // Generated Pt Vs generated y, for the cascade candidates associated with MC + Info Comb. PID
    f2dHistAsMCandCombPIDGenPtVsGenYXiMinus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYXiPlus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus(0),
    // Generated Pt Vs generated y, for the cascade candidates associated with MC
    f2dHistAsMCGenPtVsGenYXiMinus(0),
    f2dHistAsMCGenPtVsGenYXiPlus(0),
    f2dHistAsMCGenPtVsGenYOmegaMinus(0),
    f2dHistAsMCGenPtVsGenYOmegaPlus(0),
    // Generated Eta of the the cascade candidates associated with MC
    fHistAsMCGenEtaXiMinus(0),
    fHistAsMCGenEtaXiPlus(0),
    fHistAsMCGenEtaOmegaMinus(0),
    fHistAsMCGenEtaOmegaPlus(0),
    // Resolution in Pt as function of generated Pt
    f2dHistAsMCResPtXiMinus(0),		
    f2dHistAsMCResPtXiPlus(0),		
    f2dHistAsMCResPtOmegaMinus(0),
    f2dHistAsMCResPtOmegaPlus(0),	
    // Resolution in R(2D) as function of generated R
    f2dHistAsMCResRXiMinus(0),		
    f2dHistAsMCResRXiPlus(0),		
    f2dHistAsMCResROmegaMinus(0),
    f2dHistAsMCResROmegaPlus(0),
    // Resolution in phi as function of generated Pt
    f2dHistAsMCResPhiXiMinus(0),
    f2dHistAsMCResPhiXiPlus(0),
    f2dHistAsMCResPhiOmegaMinus(0),
    f2dHistAsMCResPhiOmegaPlus(0),
    // Correlation between proton (antiproton) daughter MC pt and Xi/Omega MC pt (to apply Geat/Fluka correction)
    f2dHistAsMCptProtonMCptXiMinus(0),
    f2dHistAsMCptAntiprotonMCptXiPlus(0),
    f2dHistAsMCptProtonMCptOmegaMinus(0),
    f2dHistAsMCptAntiprotonMCptOmegaPlus(0),
    // QA plots
    fHistV0toXiCosineOfPointingAngle(0),
    fHistV0CosineOfPointingAnglevsPtXi(0),
    fHistV0CosineOfPointingAnglevsPtOmega(0), 
    
    // Containers                       
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
     
        
     
//_____Non-default Constructor________________________________________________________________
AliAnalysisTaskCheckPerformanceCascadepp276::AliAnalysisTaskCheckPerformanceCascadepp276(const char *name) 
  : AliAnalysisTaskSE(name),
    fAnalysisType                  ("ESD"), 
    fESDtrackCuts                  (0),
    fPIDResponse                   (0),
    fkRerunV0CascVertexers         (0),
    fkSDDselectionOn               (kTRUE),
    fkQualityCutZprimVtxPos        (kTRUE),
    fkRejectEventPileUp            (kTRUE),
    fkQualityCutNoTPConlyPrimVtx   (kTRUE),
    fkQualityCutTPCrefit           (kTRUE),
    fkQualityCutnTPCcls            (kTRUE),
    fwithSDD                       (kTRUE),
    fMinnTPCcls                    (0),
    fkExtraSelections              (0),
    fVtxRange                      (0),
    fVtxRangeMin                   (0),
    fApplyAccCut                   (0),
    fMinPtCutOnDaughterTracks      (0),
    fEtaCutOnDaughterTracks        (0),

    // - Plots initialisation
    fListHistCascade(0),

    // - General Plots
    // Cascade multiplicity plots
    fHistCascadeMultiplicityBeforeAnySel(0),
    fHistCascadeMultiplicityAfterSDDSel(0),
    fHistCascadeMultiplicityAfterPhysicsSel(0),
    fHistCascadeMultiplicityForSelEvtNoTPCOnly(0),
    fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup(0),
    fHistCascadeMultiplicityAfterVertexCutSel(0),
    fHistnXiPlusPerEvTot(0),                  // After any event selections, in all the eta and pt range
    fHistnXiMinusPerEvTot(0),                 // After any event selections, in all the eta and pt range
    fHistnOmegaPlusPerEvTot(0),               // After any event selections, in all the eta and pt range
    fHistnOmegaMinusPerEvTot(0),              // After any event selections, in all the eta and pt range
    fHistnXiPlusPerEv(0),                     // After any event selections, in the detector acceptance and over a pt minimum
    fHistnXiMinusPerEv(0),                    // After any event selections, in the detector acceptance and over a pt minimum
    fHistnOmegaPlusPerEv(0),                  // After any event selections, in the detector acceptance and over a pt minimum
    fHistnOmegaMinusPerEv(0),                 // After any event selections, in the detector acceptance and over a pt minimum
    fHistnAssoXiMinus(0),                     // For the Reconstructed-Associated cascades 
    fHistnAssoXiPlus(0),                      // For the Reconstructed-Associated cascades 
    fHistnAssoOmegaMinus(0),                  // For the Reconstructed-Associated cascades 
    fHistnAssoOmegaPlus(0),                   // For the Reconstructed-Associated cascades 
    // Tracks multiplicity plots
    fHistTrackMultiplicityBeforeAnySel(0),
    fHistTrackMultiplicityAfterSDDSel(0),
    fHistTrackMultiplicityAfterPhysicsSel(0),
    fHistTrackMultiplicityForSelEvtNoTPCOnly(0),
    fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup(0),
    fHistTrackMultiplicityAfterVertexCutSel(0),
    // Vertex position plots (BestVertex)
    fHistPVx(0),                              // After any selections but before |Z| < 10 cm
    fHistPVy(0),                              // After any selections but before |Z| < 10 cm
    fHistPVz(0),                              // After any selections but before |Z| < 10 cm
    fHistPVxAnalysis(0),                      // After any event selections
    fHistPVyAnalysis(0),                      // After any event selections
    fHistPVzAnalysis(0),                      // After any event selections
    // - Plots before Physics Selection
    f3dHistGenPtVsGenYvsNtracksXiMinus(0),    // After the SDD event selection (For efficiency calculation)
    f3dHistGenPtVsGenctauvsYXiMinus(0),       // After the SDD event selection (For efficiency calculation)
    f3dHistGenPtVsGenYvsNtracksXiPlus(0),     // After the SDD event selection (For efficiency calculation)
    f3dHistGenPtVsGenctauvsYXiPlus(0),        // After the SDD event selection (For efficiency calculation)
    f3dHistGenPtVsGenYvsNtracksOmegaMinus(0), // After the SDD event selection (For efficiency calculation)
    f3dHistGenPtVsGenctauvsYOmegaMinus(0),    // After the SDD event selection (For efficiency calculation)
    f3dHistGenPtVsGenYvsNtracksOmegaPlus(0),  // After the SDD event selection (For efficiency calculation)
    f3dHistGenPtVsGenctauvsYOmegaPlus(0),     // After the SDD event selection (For efficiency calculation)
    // - Generated cascade plots
    // After all the event selections 
    //Xi-
    fHistEtaGenCascXiMinus(0),                // In all the eta and pt range (as they are generated)
    fHistThetaGenCascXiMinus(0),              // In all the eta and pt range (as they are generated)
    f3dHistGenPtVsGenYvsNtracksXiMinusPhysEff(0),    // 
    f3dHistGenPtVsGenctauvsYXiMinusPhysEff(0),       //
    f2dHistGenPtVsGenYFdblXiMinus(0),         // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaLambdaXiMinus(0),               // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBachXiMinus(0),                 // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaMesDghterXiMinus(0),            // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBarDghterXiMinus(0),            // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBachXiMinus(0),                    // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtMesDghterXiMinus(0),               // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBarDghterXiMinus(0),               // In the detector acceptance and over a pt minimum (Findable particle)
    //Xi+
    fHistEtaGenCascXiPlus(0),                 // In all the eta and pt range (as they are generated)
    fHistThetaGenCascXiPlus(0),               // In all the eta and pt range (as they are generated)
    f3dHistGenPtVsGenYvsNtracksXiPlusPhysEff(0),    // 
    f3dHistGenPtVsGenctauvsYXiPlusPhysEff(0),       //
    f2dHistGenPtVsGenYFdblXiPlus(0),          // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaLambdaXiPlus(0),                // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBachXiPlus(0),                  // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaMesDghterXiPlus(0),             // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBarDghterXiPlus(0),             // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBachXiPlus(0),                     // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtMesDghterXiPlus(0),                // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBarDghterXiPlus(0),                // In the detector acceptance and over a pt minimum (Findable particle)
    //Omega-
    fHistEtaGenCascOmegaMinus(0),             // In all the eta and pt range (as they are generated)
    fHistThetaGenCascOmegaMinus(0),           // In all the eta and pt range (as they are generated)
    f3dHistGenPtVsGenYvsNtracksOmegaMinusPhysEff(0),    // 
    f3dHistGenPtVsGenctauvsYOmegaMinusPhysEff(0),       //
    f2dHistGenPtVsGenYFdblOmegaMinus(0),      // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaLambdaOmegaMinus(0),            // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBachOmegaMinus(0),              // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaMesDghterOmegaMinus(0),         // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBarDghterOmegaMinus(0),         // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBachOmegaMinus(0),                 // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtMesDghterOmegaMinus(0),            // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBarDghterOmegaMinus(0),            // In the detector acceptance and over a pt minimum (Findable particle)
    //Omega+      
    fHistEtaGenCascOmegaPlus(0),              // In all the eta and pt range (as they are generated)
    fHistThetaGenCascOmegaPlus(0),            // In all the eta and pt range (as they are generated)
    f3dHistGenPtVsGenYvsNtracksOmegaPlusPhysEff(0),    // 
    f3dHistGenPtVsGenctauvsYOmegaPlusPhysEff(0),       //
    f2dHistGenPtVsGenYFdblOmegaPlus(0),       // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaLambdaOmegaPlus(0),             // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBachOmegaPlus(0),               // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaMesDghterOmegaPlus(0),          // In the detector acceptance and over a pt minimum (Findable particle)
    fHistThetaBarDghterOmegaPlus(0),          // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBachOmegaPlus(0),                  // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtMesDghterOmegaPlus(0),             // In the detector acceptance and over a pt minimum (Findable particle)
    fHistPtBarDghterOmegaPlus(0),             // In the detector acceptance and over a pt minimum (Findable particle)

    // - Associated to MC cascade plots
    fHistMassXiMinus(0),                      // For the Reconstructed-Associated cascades
    fHistMassXiPlus(0),                       // For the Reconstructed-Associated cascades
    fHistMassOmegaMinus(0),                   // For the Reconstructed-Associated cascades
    fHistMassOmegaPlus(0),                    // For the Reconstructed-Associated cascades
    // Effective mass histos with combined PID
    fHistMassWithCombPIDXiMinus(0),
    fHistMassWithCombPIDXiPlus(0),
    fHistMassWithCombPIDOmegaMinus(0),
    fHistMassWithCombPIDOmegaPlus(0),
    // PID Probability versus MC Pt(bachelor track)
    f2dHistPIDprobaKaonVsMCPtBach(0), f2dHistPIDprobaPionVsMCPtBach(0),
    // Effective mass histos with perfect MC PID on the bachelor
    fHistMassWithMcPIDXiMinus(0), fHistMassWithMcPIDXiPlus(0),
    fHistMassWithMcPIDOmegaMinus(0), fHistMassWithMcPIDOmegaPlus(0),
    // Effective mass histos for the cascade candidates associated with MC
    fHistAsMCMassXiMinus(0),
    fHistAsMCMassXiPlus(0),
    fHistAsMCMassOmegaMinus(0),
    fHistAsMCMassOmegaPlus(0),
    // Generated Pt Vs generated y, for the cascade candidates associated with MC + Info Comb. PID
    f2dHistAsMCandCombPIDGenPtVsGenYXiMinus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYXiPlus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus(0),
    // Generated Pt Vs generated y, for the cascade candidates associated with MC
    f2dHistAsMCGenPtVsGenYXiMinus(0),
    f2dHistAsMCGenPtVsGenYXiPlus(0),
    f2dHistAsMCGenPtVsGenYOmegaMinus(0),
    f2dHistAsMCGenPtVsGenYOmegaPlus(0),
    // Generated Eta of the the cascade candidates associated with MC
    fHistAsMCGenEtaXiMinus(0),
    fHistAsMCGenEtaXiPlus(0),
    fHistAsMCGenEtaOmegaMinus(0),
    fHistAsMCGenEtaOmegaPlus(0),
    // Resolution in Pt as function of generated Pt
    f2dHistAsMCResPtXiMinus(0),
    f2dHistAsMCResPtXiPlus(0),
    f2dHistAsMCResPtOmegaMinus(0),
    f2dHistAsMCResPtOmegaPlus(0),
    // Resolution in R(2D) as function of generated R
    f2dHistAsMCResRXiMinus(0),
    f2dHistAsMCResRXiPlus(0),
    f2dHistAsMCResROmegaMinus(0),
    f2dHistAsMCResROmegaPlus(0),
    // Resolution in phi as function of generated Pt
    f2dHistAsMCResPhiXiMinus(0),
    f2dHistAsMCResPhiXiPlus(0),
    f2dHistAsMCResPhiOmegaMinus(0),
    f2dHistAsMCResPhiOmegaPlus(0),
    // Correlation between proton (antiproton) daughter MC pt and Xi/Omega MC pt (to apply Geat/Fluka correction)
    f2dHistAsMCptProtonMCptXiMinus(0),
    f2dHistAsMCptAntiprotonMCptXiPlus(0),
    f2dHistAsMCptProtonMCptOmegaMinus(0),
    f2dHistAsMCptAntiprotonMCptOmegaPlus(0),
    // QA plots
    fHistV0toXiCosineOfPointingAngle(0),
    fHistV0CosineOfPointingAnglevsPtXi(0),
    fHistV0CosineOfPointingAnglevsPtOmega(0),

    // Containers                       
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
        
        // PbPb default cuts  
        fV0Sels[0] =  33.  ;     // max allowed chi2
        fV0Sels[1] =   0.1;      // min allowed impact parameter for the 1st daughter 
        fV0Sels[2] =   0.1;      // min allowed impact parameter for the 2nd daughter 
        fV0Sels[3] =   1.0 ;     // max allowed DCA between the daughter tracks       
        fV0Sels[4] =   0.998 ;   // min allowed cosine of V0's pointing angle         
        fV0Sels[5] =   0.9;      // min radius of the fiducial volume                 
        fV0Sels[6] = 100.  ;     // max radius of the fiducial volume                 
        fCascSels[0] =  33.   ;  // max allowed chi2 
        fCascSels[1] =   0.05;   // min allowed V0 impact parameter                    
        fCascSels[2] =   0.008;  // "window" around the Lambda mass                    
        fCascSels[3] =   0.03;   // min allowed bachelor's impact parameter            
        fCascSels[4] =   0.3  ;  // max allowed DCA between the V0 and the bachelor    
        fCascSels[5] =   0.999;  // min allowed cosine of the cascade pointing angle   
        fCascSels[6] =   0.9  ;  // min radius of the fiducial volume                  
        fCascSels[7] = 100.   ;  // max radius of the fiducial volume                  
              
        DefineOutput(1, TList::Class());
        DefineOutput(2, AliCFContainer::Class());
        DefineOutput(3, AliCFContainer::Class());
        DefineOutput(4, AliCFContainer::Class());
        DefineOutput(5, AliCFContainer::Class());
        DefineOutput(6, AliCFContainer::Class());
    }

    //____Destructor____
    AliAnalysisTaskCheckPerformanceCascadepp276::~AliAnalysisTaskCheckPerformanceCascadepp276()
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
void AliAnalysisTaskCheckPerformanceCascadepp276::UserCreateOutputObjects() {
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
 
 //----------------------
 // Initialize the histos
 //----------------------

 //----------------------------------
 // - Same general binning definition
 Double_t ptBinLimits[101];
 for (Int_t iptbin = 0; iptbin<101; ++iptbin) ptBinLimits[iptbin]=iptbin*0.1;
 Double_t yBinLimits[111];
 for (Int_t iybin = 0; iybin<111; ++iybin) yBinLimits[iybin]=-1.1+iybin*0.02;
 Double_t ctauBinLimits[112];
 for (Int_t ict = 0; ict<112; ++ict) ctauBinLimits[ict] = (Double_t) (ict-1.); 
 
 //------------------
 // - General plots
   // - Cascades multiplicity plots 
   if(! fHistCascadeMultiplicityBeforeAnySel) {
        fHistCascadeMultiplicityBeforeAnySel = new TH1F("fHistCascadeMultiplicityBeforeAnySel",
                        "Cascades per event (before any selections);Nbr of Cascades/Evt;Events", 50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityBeforeAnySel);
   }
   if(! fHistCascadeMultiplicityAfterSDDSel) {
        fHistCascadeMultiplicityAfterSDDSel = new TH1F("fHistCascadeMultiplicityAfterSDDSel",
                        "Cascades per event (after only the SDD selection);Nbr of Cascades/Evt;Events", 50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterSDDSel);
   }
   if(! fHistCascadeMultiplicityAfterPhysicsSel) {
        fHistCascadeMultiplicityAfterPhysicsSel = new TH1F("fHistCascadeMultiplicityAfterPhysicsSel",
                        "Cascades per event (after physics selection);Nbr of Cascades/Evt;Events", 50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterPhysicsSel);
   }
   if(! fHistCascadeMultiplicityForSelEvtNoTPCOnly) {
        fHistCascadeMultiplicityForSelEvtNoTPCOnly = new TH1F("fHistCascadeMultiplicityForSelEvtNoTPCOnly",
                        "Cascades per event (for selected events with well-established PV);Nbr of Cascades/Evt;Events", 50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityForSelEvtNoTPCOnly);
   }
   if(! fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup) {
        fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup = new TH1F("fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup",
                        "Cascades per event (for selected events with well-establisched PV and no pile-up);Nbr of Cascades/Evt;Events", 50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup);
   }
   if(! fHistCascadeMultiplicityAfterVertexCutSel) {
        fHistCascadeMultiplicityAfterVertexCutSel = new TH1F("fHistCascadeMultiplicityAfterVertexCutSel",
                                                             "Cascades per event (after vertex cut selection);Nbr of Cascades/Evt;Events", 50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterVertexCutSel);
   }
   // - Tracks multiplicity plots 
   if(! fHistTrackMultiplicityBeforeAnySel) {
        fHistTrackMultiplicityBeforeAnySel = new TH1F("fHistTrackMultiplicityBeforeAnySel",
                        "Tracks per event (before any selections);Nbr of Tracks/Evt;Events", 200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityBeforeAnySel);
   }
   if(! fHistTrackMultiplicityAfterSDDSel) {
        fHistTrackMultiplicityAfterSDDSel = new TH1F("fHistTrackMultiplicityAfterSDDSel",
                        "Tracks per event (after only the SDD selection);Nbr of Tracks/Evt;Events", 200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityAfterSDDSel);
   }
   if(! fHistTrackMultiplicityAfterPhysicsSel) {
        fHistTrackMultiplicityAfterPhysicsSel = new TH1F("fHistTrackMultiplicityAfterPhysicsSel",
                        "Tracks per event (after physics selection);Nbr of Tracks/Evt;Events", 200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityAfterPhysicsSel);
   }
   if(! fHistTrackMultiplicityForSelEvtNoTPCOnly) {
        fHistTrackMultiplicityForSelEvtNoTPCOnly = new TH1F("fHistTrackMultiplicityForSelEvtNoTPCOnly",
                        "Tracks per event (for selected events with well-established PV);Nbr of Tracks/Evt;Events", 200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityForSelEvtNoTPCOnly);
   }
   if(! fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup) {
        fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup = new TH1F("fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup",
                        "Tracks per event (for selected events with well-establisched PV and no pile-up);Nbr of Tracks/Evt;Events", 200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup);
   }
   if(! fHistTrackMultiplicityAfterVertexCutSel) {
        fHistTrackMultiplicityAfterVertexCutSel = new TH1F("fHistTrackMultiplicityAfterVertexCutSel",
                                                           "Tracks per event (after vertex cut selection);Nbr of Tracks/Evt;Events", 200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityAfterVertexCutSel);
   }
   // - Vertex position plots
   if(! fHistPVx ){
        fHistPVx = new TH1F("fHistPVx", "Best PV position in x; x (cm); Events", 2000, -0.5, 0.5);
        fListHistCascade->Add(fHistPVx);
   }
   if(! fHistPVy ){
        fHistPVy = new TH1F("fHistPVy", "Best PV position in y; y (cm); Events", 2000, -0.5, 0.5);
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
        fHistPVyAnalysis = new TH1F("fHistPVyAnalysis", "Best PV position in y (after events selections); y (cm); Events" , 2000, -0.5, 0.5);
        fListHistCascade->Add(fHistPVyAnalysis);
   }
   if(! fHistPVzAnalysis ){
        fHistPVzAnalysis = new TH1F("fHistPVzAnalysis", "Best PV position in z (after events selections); z (cm); Events", 400, -20, 20);
        fListHistCascade->Add(fHistPVzAnalysis);
   }

 //--------------------------
 // - Generated cascade plots
   // - Generated Cascade multiplicity distributions (for singol cascade)
   fHistnXiPlusPerEvTot = new TH1F("fHistnXiPlusPerEvTot", "", 25, 0, 25);
   fListHistCascade->Add(fHistnXiPlusPerEvTot);
   fHistnXiMinusPerEvTot = new TH1F("fHistnXiMinusPerEvTot", "", 25, 0, 25);
   fListHistCascade->Add(fHistnXiMinusPerEvTot);
   fHistnOmegaPlusPerEvTot = new TH1F("fHistnOmegaPlusPerEvTot", "", 25, 0, 25);
   fListHistCascade->Add(fHistnOmegaPlusPerEvTot);
   fHistnOmegaMinusPerEvTot = new TH1F("fHistnOmegaMinusPerEvTot", "", 25, 0, 25);
   fListHistCascade->Add(fHistnOmegaMinusPerEvTot);   
   fHistnXiPlusPerEv = new TH1F("fHistnXiPlusPerEv", "", 25, 0, 25);
   fListHistCascade->Add(fHistnXiPlusPerEv);
   fHistnXiMinusPerEv = new TH1F("fHistnXiMinusPerEv", "", 25, 0, 25);
   fListHistCascade->Add(fHistnXiMinusPerEv);
   fHistnOmegaPlusPerEv = new TH1F("fHistnOmegaPlusPerEv", "", 25, 0, 25);
   fListHistCascade->Add(fHistnOmegaPlusPerEv);
   fHistnOmegaMinusPerEv = new TH1F("fHistnOmegaMinusPerEv", "", 25, 0, 25);
   fListHistCascade->Add(fHistnOmegaMinusPerEv);
   // - Xi- 
   // - Pseudo-Rapidity distribution
   if (!fHistEtaGenCascXiMinus) {
     fHistEtaGenCascXiMinus = new TH1F("fHistEtaGenCascXiMinus", "#eta of any gen. #Xi^{-}; #eta; Number of Casc", 200, -10, 10);
     fListHistCascade->Add(fHistEtaGenCascXiMinus);
   }
   if (!f3dHistGenPtVsGenYvsNtracksXiMinus) {
     f3dHistGenPtVsGenYvsNtracksXiMinus = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinus", "MC P_{t} Vs MC Y of Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinus);
   }
   if (!f3dHistGenPtVsGenctauvsYXiMinus) {
      f3dHistGenPtVsGenctauvsYXiMinus = new TH3D("f3dHistGenPtVsGenctauvsYXiMinus", "MC P_{t} Vs MC ctau Vs Y of Gen #Xi^{-}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
      fListHistCascade->Add(f3dHistGenPtVsGenctauvsYXiMinus);
   }
   if (!f3dHistGenPtVsGenYvsNtracksXiMinusPhysEff) {
     f3dHistGenPtVsGenYvsNtracksXiMinusPhysEff = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinusPhysEff", "MC P_{t} Vs MC Y of Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinusPhysEff);
   }
   if (!f3dHistGenPtVsGenctauvsYXiMinusPhysEff) {
      f3dHistGenPtVsGenctauvsYXiMinusPhysEff = new TH3D("f3dHistGenPtVsGenctauvsYXiMinusPhysEff", "MC P_{t} Vs MC ctau Vs Y of Gen #Xi^{-}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
      fListHistCascade->Add(f3dHistGenPtVsGenctauvsYXiMinusPhysEff);
   }
   // - Info at the generation level of multi-strange particle
   if (!fHistThetaGenCascXiMinus) {
      fHistThetaGenCascXiMinus = new TH1F("fHistThetaGenCascXiMinus", "#theta of gen. #Xi^{-}; #theta; Number of Casc.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaGenCascXiMinus);
   }
   if (!f2dHistGenPtVsGenYFdblXiMinus) {
      f2dHistGenPtVsGenYFdblXiMinus = new TH2D("f2dHistGenPtVsGenYFdblXiMinus", "MC P_{t} Vs MC Y of findable Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistGenPtVsGenYFdblXiMinus);
   }
   // - Theta distribution the daughters (control plots)
   if (!fHistThetaLambdaXiMinus) {
      fHistThetaLambdaXiMinus = new TH1F("fHistThetaLambdaXiMinus", "#theta of gen. #Lambda (Xi dghter); #theta_{#Lambda}; Number of #Lambda^0", 200, -10, 190);
      fListHistCascade->Add(fHistThetaLambdaXiMinus);
   }
   if (!fHistThetaBachXiMinus) {
      fHistThetaBachXiMinus = new TH1F("fHistThetaBachXiMinus", "#theta of gen. Bach.; #theta_{Bach}; Number of Bach.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaBachXiMinus);
   }
   if (!fHistThetaMesDghterXiMinus) {
      fHistThetaMesDghterXiMinus = new TH1F("fHistThetaMesDghterXiMinus", "#theta of gen. Meson #Lambda dghter; #theta_{MesDght}; Number of Mes.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaMesDghterXiMinus);
   }
   if (!fHistThetaBarDghterXiMinus) {
      fHistThetaBarDghterXiMinus = new TH1F("fHistThetaBarDghterXiMinus", "#theta of gen. Baryon #Lambda dghter; #theta_{BarDght}; Number of Bar.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaBarDghterXiMinus);
   }
   // - Pt distribution (control plots)
   if (!fHistPtBachXiMinus) {
      fHistPtBachXiMinus = new TH1F("fHistPtBachXiMinus", "p_{t} of gen. Bach.; pt_{Bach}; Number of Bach.", 200, 0, 10);
      fListHistCascade->Add(fHistPtBachXiMinus);
   }
   if (!fHistPtMesDghterXiMinus) {
      fHistPtMesDghterXiMinus = new TH1F("fHistPtMesDghterXiMinus", "p_{t} of gen. Meson #Lambda dghter; pt_{MesDght}; Number of Mes.", 200, 0, 10);
      fListHistCascade->Add(fHistPtMesDghterXiMinus);
   }
   if (!fHistPtBarDghterXiMinus) {
      fHistPtBarDghterXiMinus = new TH1F("fHistPtBarDghterXiMinus", "p_{t} of gen. Baryon #Lambda dghter; pt_{BarDght}; Number of Bar.", 200, 0, 10);
      fListHistCascade->Add(fHistPtBarDghterXiMinus);
   }
   // - Xi+ 
   // - Pseudo-Rapidity distribution
   if (!fHistEtaGenCascXiPlus) {
      fHistEtaGenCascXiPlus = new TH1F("fHistEtaGenCascXiPlus", "#eta of any gen. #Xi^{+}; #eta; Number of Casc", 200, -10, 10);
      fListHistCascade->Add(fHistEtaGenCascXiPlus);
   }
   if (!f3dHistGenPtVsGenYvsNtracksXiPlus) {
      f3dHistGenPtVsGenYvsNtracksXiPlus = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlus", "MC P_{t} Vs MC Y of Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
      fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlus);
   }
   if (!f3dHistGenPtVsGenctauvsYXiPlus) {
      f3dHistGenPtVsGenctauvsYXiPlus = new TH3D("f3dHistGenPtVsGenctauvsYXiPlus", "MC P_{t} Vs MC ctau Vs Yof Gen #Xi^{+}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
      fListHistCascade->Add(f3dHistGenPtVsGenctauvsYXiPlus);
   }
   if (!f3dHistGenPtVsGenYvsNtracksXiPlusPhysEff) {
      f3dHistGenPtVsGenYvsNtracksXiPlusPhysEff = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlusPhysEff", "MC P_{t} Vs MC Y of Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
      fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlusPhysEff);
   }
   if (!f3dHistGenPtVsGenctauvsYXiPlusPhysEff) {
      f3dHistGenPtVsGenctauvsYXiPlusPhysEff = new TH3D("f3dHistGenPtVsGenctauvsYXiPlusPhysEff", "MC P_{t} Vs MC ctau Vs Yof Gen #Xi^{+}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
      fListHistCascade->Add(f3dHistGenPtVsGenctauvsYXiPlusPhysEff);
   }
   // - Info at the generation level of multi-strange particle
   if (!fHistThetaGenCascXiPlus) {
      fHistThetaGenCascXiPlus = new TH1F("fHistThetaGenCascXiPlus", "#theta of gen. #Xi^{+}; #theta; Number of Casc.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaGenCascXiPlus);
   }
   if (!f2dHistGenPtVsGenYFdblXiPlus) {
      f2dHistGenPtVsGenYFdblXiPlus = new TH2D("f2dHistGenPtVsGenYFdblXiPlus", "MC P_{t} Vs MC Y of findable Gen #Xi^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistGenPtVsGenYFdblXiPlus);
   }
   // - Theta distribution the daughters (control plots)
   if (!fHistThetaLambdaXiPlus) {
      fHistThetaLambdaXiPlus = new TH1F("fHistThetaLambdaXiPlus", "#theta of gen. #Lambda (Xi dghter); #theta_{#Lambda}; Number of #Lambda", 200, -10, 190);
      fListHistCascade->Add(fHistThetaLambdaXiPlus);
   }
   if (!fHistThetaBachXiPlus) {
      fHistThetaBachXiPlus = new TH1F("fHistThetaBachXiPlus", "#theta of gen. Bach.; #theta_{Bach}; Number of Bach.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaBachXiPlus);
   }
   if (!fHistThetaMesDghterXiPlus) {
      fHistThetaMesDghterXiPlus = new TH1F("fHistThetaMesDghterXiPlus", "#theta of gen. Meson #Lambda dghter; #theta_{MesDght}; Number of Mes.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaMesDghterXiPlus);
   }
   if (!fHistThetaBarDghterXiPlus) {
      fHistThetaBarDghterXiPlus = new TH1F("fHistThetaBarDghterXiPlus", "#theta of gen. Baryon #Lambda dghter; #theta_{BarDght}; Number of Bar.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaBarDghterXiPlus);
   }
   // - Pt distribution (control plots)
   if (!fHistPtBachXiPlus) {
      fHistPtBachXiPlus = new TH1F("fHistPtBachXiPlus", "p_{t} of gen. Bach.; pt_{Bach}; Number of Bach.", 200, 0, 10);
      fListHistCascade->Add(fHistPtBachXiPlus);
   }
   if (!fHistPtMesDghterXiPlus) {
      fHistPtMesDghterXiPlus = new TH1F("fHistPtMesDghterXiPlus", "p_{t} of gen. Meson #Lambda dghter; pt_{MesDght}; Number of Mes.", 200, 0, 10);
      fListHistCascade->Add(fHistPtMesDghterXiPlus);
   }
   if (!fHistPtBarDghterXiPlus) {
      fHistPtBarDghterXiPlus = new TH1F("fHistPtBarDghterXiPlus", "p_{t} of gen. Baryon #Lambda dghter); pt_{BarDght}; Number of Bar.", 200, 0, 10);
      fListHistCascade->Add(fHistPtBarDghterXiPlus);
   }
   // - Omega- 
   // - Pseudo-Rapidity distribution
   if (!fHistEtaGenCascOmegaMinus) {
      fHistEtaGenCascOmegaMinus = new TH1F("fHistEtaGenCascOmegaMinus", "#eta of any gen. #Omega^{-}; #eta; Number of Casc", 200, -10, 10);
      fListHistCascade->Add(fHistEtaGenCascOmegaMinus);
   }
   if (!f3dHistGenPtVsGenYvsNtracksOmegaMinus) {
      f3dHistGenPtVsGenYvsNtracksOmegaMinus = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinus", "MC P_{t} Vs MC Y of Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
      fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinus);
   }
   if (!f3dHistGenPtVsGenctauvsYOmegaMinus) {
      f3dHistGenPtVsGenctauvsYOmegaMinus = new TH3D("f3dHistGenPtVsGenctauvsYOmegaMinus", "MC P_{t} Vs MC ctau Vs Y of Gen #Omega^{-} ", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
      fListHistCascade->Add(f3dHistGenPtVsGenctauvsYOmegaMinus);
   }
   if (!f3dHistGenPtVsGenYvsNtracksOmegaMinusPhysEff) {
      f3dHistGenPtVsGenYvsNtracksOmegaMinusPhysEff = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinusPhysEff", "MC P_{t} Vs MC Y of Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
      fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinusPhysEff);
   }
   if (!f3dHistGenPtVsGenctauvsYOmegaMinusPhysEff) {
      f3dHistGenPtVsGenctauvsYOmegaMinusPhysEff = new TH3D("f3dHistGenPtVsGenctauvsYOmegaMinusPhysEff", "MC P_{t} Vs MC ctau Vs Y of Gen #Omega^{-}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
      fListHistCascade->Add(f3dHistGenPtVsGenctauvsYOmegaMinusPhysEff);
   }
   // - Info at the generation level of multi-strange particle
   if (!fHistThetaGenCascOmegaMinus) {
      fHistThetaGenCascOmegaMinus = new TH1F("fHistThetaGenCascOmegaMinus", "#theta of gen. #Omega^{-}; #theta; Number of Casc.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaGenCascOmegaMinus);
   }
   if (!f2dHistGenPtVsGenYFdblOmegaMinus) {
      f2dHistGenPtVsGenYFdblOmegaMinus = new TH2D("f2dHistGenPtVsGenYFdblOmegaMinus", "MC P_{t} Vs MC Y of findable Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistGenPtVsGenYFdblOmegaMinus);
   }
   // - Theta distribution the daughters (control plots)
   if (!fHistThetaLambdaOmegaMinus) {
      fHistThetaLambdaOmegaMinus = new TH1F("fHistThetaLambdaOmegaMinus", "#theta of gen. #Lambda (Omega dghter); #theta_{#Lambda}; Number of #Lambda", 200, -10, 190);
      fListHistCascade->Add(fHistThetaLambdaOmegaMinus);
   }
   if (!fHistThetaBachOmegaMinus) {
      fHistThetaBachOmegaMinus = new TH1F("fHistThetaBachOmegaMinus", "#theta of gen. Bach.;#theta_{Bach};Number of Bach.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaBachOmegaMinus);
   }
   if (!fHistThetaMesDghterOmegaMinus) {
      fHistThetaMesDghterOmegaMinus = new TH1F("fHistThetaMesDghterOmegaMinus", "#theta of gen. Meson #Lambda dghter; #theta_{MesDght}; Number of Mes.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaMesDghterOmegaMinus);
   }
   if (!fHistThetaBarDghterOmegaMinus) {
      fHistThetaBarDghterOmegaMinus = new TH1F("fHistThetaBarDghterOmegaMinus", "#theta of gen. Baryon #Lambda dghter; #theta_{BarDght}; Number of Bar.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaBarDghterOmegaMinus);
   }
   // - Pt distribution (control plots)
   if (!fHistPtBachOmegaMinus) {
      fHistPtBachOmegaMinus = new TH1F("fHistPtBachOmegaMinus", "p_{t} of gen. Bach.; pt_{Bach}; Number of Bach.", 200, 0, 10);
      fListHistCascade->Add(fHistPtBachOmegaMinus);
   }
   if (!fHistPtMesDghterOmegaMinus) {
      fHistPtMesDghterOmegaMinus = new TH1F("fHistPtMesDghterOmegaMinus", "p_{t} of gen. Meson #Lambda dghter); pt_{MesDght}; Number of Mes.", 200, 0, 10);
      fListHistCascade->Add(fHistPtMesDghterOmegaMinus);
   }
   if (!fHistPtBarDghterOmegaMinus) {
      fHistPtBarDghterOmegaMinus = new TH1F("fHistPtBarDghterOmegaMinus", "p_{t} of gen. Baryon #Lambda dghter); pt_{BarDght}; Number of Bar.", 200, 0, 10);
      fListHistCascade->Add(fHistPtBarDghterOmegaMinus);
   }
   // - Omega+ 
   // - Pseudo-Rapidity distribution
   if (!fHistEtaGenCascOmegaPlus) {
      fHistEtaGenCascOmegaPlus = new TH1F("fHistEtaGenCascOmegaPlus", "#eta of any gen. #Omega^{+}; #eta; Number of Casc", 200, -10, 10);
      fListHistCascade->Add(fHistEtaGenCascOmegaPlus);
   }
   if (!f3dHistGenPtVsGenYvsNtracksOmegaPlus) {
      f3dHistGenPtVsGenYvsNtracksOmegaPlus = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlus", "MC P_{t} Vs MC Y of Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
      fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlus);
   }
   if (!f3dHistGenPtVsGenctauvsYOmegaPlus) {
      f3dHistGenPtVsGenctauvsYOmegaPlus = new TH3D("f3dHistGenPtVsGenctauvsYOmegaPlus", "MC P_{t} Vs MC ctau Vs Y of Gen #Omega^{+} ", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
      fListHistCascade->Add(f3dHistGenPtVsGenctauvsYOmegaPlus);
   }
   if (!f3dHistGenPtVsGenYvsNtracksOmegaPlusPhysEff) {
      f3dHistGenPtVsGenYvsNtracksOmegaPlusPhysEff = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlusPhysEff", "MC P_{t} Vs MC Y of Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 200, 0., 200.);
      fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlusPhysEff);
   }
   if (!f3dHistGenPtVsGenctauvsYOmegaPlusPhysEff) {
      f3dHistGenPtVsGenctauvsYOmegaPlusPhysEff = new TH3D("f3dHistGenPtVsGenctauvsYOmegaPlusPhysEff", "MC P_{t} Vs MC ctau Vs Y of Gen #Omega^{+}", 100, ptBinLimits, 111, ctauBinLimits, 110, yBinLimits);
      fListHistCascade->Add(f3dHistGenPtVsGenctauvsYOmegaPlusPhysEff);
   }
   // - Info at the generation level of multi-strange particle
   if (!fHistThetaGenCascOmegaPlus) {
      fHistThetaGenCascOmegaPlus = new TH1F("fHistThetaGenCascOmegaPlus", "#theta of gen. #Omega^{+}; #theta; Number of Casc.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaGenCascOmegaPlus);
   }
   if (!f2dHistGenPtVsGenYFdblOmegaPlus) {
      f2dHistGenPtVsGenYFdblOmegaPlus = new TH2D("f2dHistGenPtVsGenYFdblOmegaPlus", "MC P_{t} Vs MC Y of findable Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistGenPtVsGenYFdblOmegaPlus);
   }
   // - Info at the generation level of multi-strange particle
   if (!fHistThetaGenCascOmegaPlus) {
      fHistThetaGenCascOmegaPlus = new TH1F("fHistThetaGenCascOmegaPlus", "#theta of gen. #Omega^{+}; #theta; Number of Casc.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaGenCascOmegaPlus);
   }
   if (!f2dHistGenPtVsGenYFdblOmegaPlus) {
      f2dHistGenPtVsGenYFdblOmegaPlus = new TH2D("f2dHistGenPtVsGenYFdblOmegaPlus", "MC P_{t} Vs MC Y of findable Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistGenPtVsGenYFdblOmegaPlus);
   }
   // - Theta distribution the daughters (control plots)
   if (!fHistThetaLambdaOmegaPlus) {
      fHistThetaLambdaOmegaPlus = new TH1F("fHistThetaLambdaOmegaPlus", "#theta of gen. #Lambda (Omega dghter); #theta_{#Lambda}; Number of #Lambda", 200, -10, 190);
      fListHistCascade->Add(fHistThetaLambdaOmegaPlus);
   }
   if (!fHistThetaBachOmegaPlus) {
      fHistThetaBachOmegaPlus = new TH1F("fHistThetaBachOmegaPlus", "#theta of gen. Bach.; #theta_{Bach}; Number of Bach.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaBachOmegaPlus);
   }
   if (!fHistThetaMesDghterOmegaPlus) {
      fHistThetaMesDghterOmegaPlus = new TH1F("fHistThetaMesDghterOmegaPlus", "#theta of gen. Meson #Lambda dghter; #theta_{MesDght}; Number of Mes.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaMesDghterOmegaPlus);
   }
   if (!fHistThetaBarDghterOmegaPlus) {
      fHistThetaBarDghterOmegaPlus = new TH1F("fHistThetaBarDghterOmegaPlus", "#theta of gen. Baryon #Lambda dghter; #theta_{BarDght}; Number of Bar.", 200, -10, 190);
      fListHistCascade->Add(fHistThetaBarDghterOmegaPlus);
   }
   // - Pt distribution (control plots)
   if (!fHistPtBachOmegaPlus) {
      fHistPtBachOmegaPlus = new TH1F("fHistPtBachOmegaPlus", "p_{t} of gen. Bach.; pt_{Bach}; Number of Bach.", 200, 0, 10);
      fListHistCascade->Add(fHistPtBachOmegaPlus);
   }
   if (!fHistPtMesDghterOmegaPlus) {
      fHistPtMesDghterOmegaPlus = new TH1F("fHistPtMesDghterOmegaPlus", "p_{t} of gen. Meson #Lambda dghter; pt_{MesDght}; Number of Mes.", 200, 0, 10);
      fListHistCascade->Add(fHistPtMesDghterOmegaPlus);
   }
   if (!fHistPtBarDghterOmegaPlus) {
      fHistPtBarDghterOmegaPlus = new TH1F("fHistPtBarDghterOmegaPlus", "p_{t} of gen. Baryon #Lambda dghter); pt_{BarDght}; Number of Bar.", 200, 0, 10);
      fListHistCascade->Add(fHistPtBarDghterOmegaPlus);
   }
 
 //-------------------------------------------------------------------------
 // - Any reconstructed cascades + reconstructed cascades associated with MC
  
   // - Multiplicity cascde plots
   fHistnAssoXiMinus= new TH1F("fHistnAssoXiMinus", "", 25, 0, 25);
   fListHistCascade->Add(fHistnAssoXiMinus);
   fHistnAssoXiPlus= new TH1F("fHistnAssoXiPlus", "", 25, 0, 25);
   fListHistCascade->Add(fHistnAssoXiPlus);
   fHistnAssoOmegaMinus= new TH1F("fHistnAssoOmegaMinus", "", 25, 0, 25);
   fListHistCascade->Add(fHistnAssoOmegaMinus);
   fHistnAssoOmegaPlus= new TH1F("fHistnAssoOmegaPlus", "", 25, 0, 25);
   fListHistCascade->Add(fHistnAssoOmegaPlus);
   // - Effective mass histos for cascades candidates. 
   if (! fHistMassXiMinus) {
	  fHistMassXiMinus = new TH1F("fHistMassXiMinus","#Xi^{-} candidates; M( #Lambda , #pi^{-} ) (GeV/c^{2}); Counts", 400, 1.2, 2.0);
 	  fListHistCascade->Add(fHistMassXiMinus);
   }
   if (! fHistMassXiPlus) {
	  fHistMassXiPlus = new TH1F("fHistMassXiPlus","#Xi^{+} candidates; M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2}); Counts", 400, 1.2, 2.0);
	  fListHistCascade->Add(fHistMassXiPlus);
   }
   if (! fHistMassOmegaMinus) {
	  fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus","#Omega^{-} candidates; M( #Lambda , K^{-} ) (GeV/c^{2}); Counts", 500, 1.5, 2.5);
	  fListHistCascade->Add(fHistMassOmegaMinus);
   } 
   if (! fHistMassOmegaPlus) {
	  fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus","#Omega^{+} candidates; M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2}); Counts", 500, 1.5, 2.5);
	  fListHistCascade->Add(fHistMassOmegaPlus);
   }
   // - Effective mass histos with combined PID
   if (! fHistMassWithCombPIDXiMinus) {
      fHistMassWithCombPIDXiMinus = new TH1F("fHistMassWithCombPIDXiMinus","#Xi^{-} candidates, with Bach. comb. PID; M( #Lambda , #pi^{-} ) (GeV/c^{2}); Counts", 400, 1.2, 2.0);
      fListHistCascade->Add(fHistMassWithCombPIDXiMinus);
   }
   if (! fHistMassWithCombPIDXiPlus) {
      fHistMassWithCombPIDXiPlus = new TH1F("fHistMassWithCombPIDXiPlus","#Xi^{+} candidates, with Bach. comb. PID; M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2}); Counts", 400, 1.2, 2.0);
      fListHistCascade->Add(fHistMassWithCombPIDXiPlus);
   }
   if (! fHistMassWithCombPIDOmegaMinus) {
      fHistMassWithCombPIDOmegaMinus = new TH1F("fHistMassWithCombPIDOmegaMinus","#Omega^{-} candidates, with Bach. comb. PID; M( #Lambda , K^{-} ) (GeV/c^{2}); Counts", 500, 1.5, 2.5);
      fListHistCascade->Add(fHistMassWithCombPIDOmegaMinus);
   }
   if (! fHistMassWithCombPIDOmegaPlus) {
      fHistMassWithCombPIDOmegaPlus = new TH1F("fHistMassWithCombPIDOmegaPlus","#Omega^{+} candidates, with Bach. comb. PID; M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2}); Counts", 500, 1.5, 2.5);
      fListHistCascade->Add(fHistMassWithCombPIDOmegaPlus);
   }
   // - PID Probability versus MC Pt(bachelor track)
   if (! f2dHistPIDprobaKaonVsMCPtBach ){
      f2dHistPIDprobaKaonVsMCPtBach  = new TH2F("f2dHistPIDprobaKaonVsMCPtBach", "Comb. PID proba to be K^{#pm} Vs MC Bach. Pt; Pt_{MC}(Bach.) (GeV/c); Comb. PID Proba (Bach. = K^{#pm})", 100, 0.0, 5.0, 110, 0.0, 1.10);
      fListHistCascade->Add(f2dHistPIDprobaKaonVsMCPtBach);
   }
   if(! f2dHistPIDprobaPionVsMCPtBach ){
      f2dHistPIDprobaPionVsMCPtBach  = new TH2F("f2dHistPIDprobaPionVsMCPtBach", "Comb. PID proba to be #pi^{#pm} Vs MC Bach. Pt; Pt_{MC}(Bach.) (GeV/c); Comb. PID Proba (Bach. = #pi^{#pm})", 100, 0.0, 5.0, 110, 0.0, 1.10);
      fListHistCascade->Add(f2dHistPIDprobaPionVsMCPtBach);
   }
   // - Effective mass histos with perfect MC PID on the bachelor
   if (! fHistMassWithMcPIDXiMinus) {
      fHistMassWithMcPIDXiMinus = new TH1F("fHistMassWithMcPIDXiMinus", "#Xi^{-} candidates, with Bach. MC PID; M( #Lambda , #pi^{-} ) (GeV/c^{2}); Counts", 400, 1.2, 2.0);
      fListHistCascade->Add(fHistMassWithMcPIDXiMinus);
   }
   if (! fHistMassWithMcPIDXiPlus) {
      fHistMassWithMcPIDXiPlus = new TH1F("fHistMassWithMcPIDXiPlus", "#Xi^{+} candidates, with Bach. MC PID; M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2}); Counts", 400, 1.2, 2.0);
      fListHistCascade->Add(fHistMassWithMcPIDXiPlus);
   }
   if (! fHistMassWithMcPIDOmegaMinus) {
      fHistMassWithMcPIDOmegaMinus = new TH1F("fHistMassWithMcPIDOmegaMinus", "#Omega^{-} candidates, with Bach. MC PID; M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500, 1.5, 2.5);
      fListHistCascade->Add(fHistMassWithMcPIDOmegaMinus);
   }
   if (! fHistMassWithMcPIDOmegaPlus) {
      fHistMassWithMcPIDOmegaPlus = new TH1F("fHistMassWithMcPIDOmegaPlus", "#Omega^{+} candidates, with Bach. MC PID; M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2}); Counts", 500, 1.5, 2.5);
      fListHistCascade->Add(fHistMassWithMcPIDOmegaPlus);
   }
   // - Effective mass histos for cascades candidates ASSOCIATED with MC.
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
   // -  Generated Pt Vs generated Y of the cascade candidates associated with MC + having the proper maximum proba of combined PID for the bachelor
   if (!f2dHistAsMCandCombPIDGenPtVsGenYXiMinus) {
      f2dHistAsMCandCombPIDGenPtVsGenYXiMinus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYXiMinus", "MC P_{t} Vs MC Y of #Xi^{-} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYXiMinus);
   }
   if (!f2dHistAsMCandCombPIDGenPtVsGenYXiPlus) {
      f2dHistAsMCandCombPIDGenPtVsGenYXiPlus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYXiPlus", "MC P_{t} Vs MC Y of #Xi^{+} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYXiPlus);
   } 
   if (!f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus) {
      f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus", "MC P_{t} Vs MC Y of #Omega^{-} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus);
   }
   if (!f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus) {
      f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus", "MC P_{t} Vs MC Y of #Omega^{+} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus);
   }
   // - Generated Pt Vs Generated Y, for the cascade candidates associated with MC
   if (!f2dHistAsMCGenPtVsGenYXiMinus) {
      f2dHistAsMCGenPtVsGenYXiMinus = new TH2F("f2dHistAsMCGenPtVsGenYXiMinus", "MC P_{t} Vs MC Y of gen. #Xi^{-} (associated); Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistAsMCGenPtVsGenYXiMinus );
   }
   if (!f2dHistAsMCGenPtVsGenYXiPlus) {
      f2dHistAsMCGenPtVsGenYXiPlus = new TH2F("f2dHistAsMCGenPtVsGenYXiPlus", "MC P_{t} Vs MC Y of gen. #Xi^{+} (associated); Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistAsMCGenPtVsGenYXiPlus );
   }
   if (!f2dHistAsMCGenPtVsGenYOmegaMinus) {
      f2dHistAsMCGenPtVsGenYOmegaMinus = new TH2F("f2dHistAsMCGenPtVsGenYOmegaMinus", "MC P_{t} Vs MC Y of gen. #Omega^{-} (associated); Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistAsMCGenPtVsGenYOmegaMinus );
   }
   if (!f2dHistAsMCGenPtVsGenYOmegaPlus) {
      f2dHistAsMCGenPtVsGenYOmegaPlus = new TH2F("f2dHistAsMCGenPtVsGenYOmegaPlus", "MC P_{t} Vs MC Y of gen. #Omega^{+} (associated); Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
      fListHistCascade->Add(f2dHistAsMCGenPtVsGenYOmegaPlus );
   } 
   // - Generated Eta of the the cascade candidates associated with MC
   if (!fHistAsMCGenEtaXiMinus) {
      fHistAsMCGenEtaXiMinus = new TH1F("fHistAsMCGenEtaXiMinus", "#eta of gen. #Xi^{-} (associated); #eta; Count", 100, -5, 5);
      fListHistCascade->Add( fHistAsMCGenEtaXiMinus );
   }
   if (!fHistAsMCGenEtaXiPlus) {
      fHistAsMCGenEtaXiPlus = new TH1F("fHistAsMCGenEtaXiPlus", "#eta of gen. #Xi^{+} (associated); #eta; Count", 100, -5, 5);
      fListHistCascade->Add( fHistAsMCGenEtaXiPlus );
   }
   if (!fHistAsMCGenEtaOmegaMinus) {
      fHistAsMCGenEtaOmegaMinus = new TH1F("fHistAsMCGenEtaOmegaMinus", "#eta of gen. #Omega^{-} (associated);#eta;Number of Casc", 100, -5, 5);
      fListHistCascade->Add( fHistAsMCGenEtaOmegaMinus );
   }
   if (!fHistAsMCGenEtaOmegaPlus) {
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
   if (!f2dHistAsMCptProtonMCptXiMinus) {
      f2dHistAsMCptProtonMCptXiMinus = new TH2F("f2dHistAsMCptProtonMCptXiMinus", "Proton MC pt vs Xi- MC pt", 100, 0., 10., 100, 0., 10.); 
      fListHistCascade->Add(f2dHistAsMCptProtonMCptXiMinus);
   }
   if (!f2dHistAsMCptAntiprotonMCptXiPlus) {
      f2dHistAsMCptAntiprotonMCptXiPlus = new TH2F("f2dHistAsMCptAntiprotonMCptXiPlus", "Antiproton MC pt vs Xi+ MC pt", 100, 0., 10., 100, 0., 10.);
      fListHistCascade->Add(f2dHistAsMCptAntiprotonMCptXiPlus);
   }
   if (!f2dHistAsMCptProtonMCptOmegaMinus) {
      f2dHistAsMCptProtonMCptOmegaMinus = new TH2F("f2dHistAsMCptProtonMCptOmegaMinus", "Proton MC pt vs Omega- MC pt", 100, 0., 10., 100, 0., 10.);
      fListHistCascade->Add(f2dHistAsMCptProtonMCptOmegaMinus);
   }
   if (!f2dHistAsMCptAntiprotonMCptOmegaPlus) {
      f2dHistAsMCptAntiprotonMCptOmegaPlus = new TH2F("f2dHistAsMCptAntiprotonMCptOmegaPlus", "Antiproton MC pt vs Omega+ MC pt", 100, 0., 10., 100, 0., 10.);
      fListHistCascade->Add(f2dHistAsMCptAntiprotonMCptOmegaPlus);
   }
   // - Cosine of Pointing angle
   if (! fHistV0toXiCosineOfPointingAngle) {
      fHistV0toXiCosineOfPointingAngle = new TH1F("fHistV0toXiCosineOfPointingAngle", "Cos. of V0 Ptng Angl / Xi vtx ; Cos(V0 Point. Angl / Xi vtx); Counts", 200, 0.95, 1.0001);
      fListHistCascade->Add(fHistV0toXiCosineOfPointingAngle);
   }
   if (! fHistV0CosineOfPointingAnglevsPtXi) {
      fHistV0CosineOfPointingAnglevsPtXi = new TH2F("fHistV0CosineOfPointingAnglevsPtXi", "Cos. of V0 Ptng Angl vs cascade Pt; Cos(V0 Point. Angl); Counts", 100, 0., 10., 200, 0.95, 1.0001);
      fListHistCascade->Add(fHistV0CosineOfPointingAnglevsPtXi);
   }
   if (! fHistV0CosineOfPointingAnglevsPtOmega) {
      fHistV0CosineOfPointingAnglevsPtOmega = new TH2F("fHistV0CosineOfPointingAnglevsPtOmega", "Cos. of V0 Ptng Angl vs cascade Pt; Cos(V0 Point. Angl); Counts", 100, 0., 10., 200, 0.95, 1.0001);
      fListHistCascade->Add(fHistV0CosineOfPointingAnglevsPtOmega);
   }

  //--------------
  // - CFContainer
  // PID container Xi-
  if(! fCFContCascadePIDAsXiMinus)  {
     const Int_t  lNbSteps      =  7;
     const Int_t  lNbVariables  =  3;
       //Array for the number of bins in each dimension:
     Int_t lNbBinsPerVar[3] = {0};
     lNbBinsPerVar[0] = 100;
     lNbBinsPerVar[1] = 800;
     lNbBinsPerVar[2] = 22;
     fCFContCascadePIDAsXiMinus = new AliCFContainer(Form("fCFContCascadePIDAsXiMinus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
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
     fListHistCascade->Add(fCFContCascadePIDAsXiMinus);  
  }
  // PID container Xi+
  if(! fCFContCascadePIDAsXiPlus)  {
     const Int_t  lNbSteps      =  7;
     const Int_t  lNbVariables  =  3;
       //Array for the number of bins in each dimension :
     Int_t lNbBinsPerVar[3] = {0};
     lNbBinsPerVar[0] = 100;
     lNbBinsPerVar[1] = 800;
     lNbBinsPerVar[2] = 22;
     fCFContCascadePIDAsXiPlus = new AliCFContainer(Form("fCFContCascadePIDAsXiPlus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
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
     fListHistCascade->Add(fCFContCascadePIDAsXiPlus);
  }
  // PID container Omega-
  if(! fCFContCascadePIDAsOmegaMinus)  {
     const Int_t  lNbSteps      =  7;
     const Int_t  lNbVariables  =  3;
       //Array for the number of bins in each dimension :
     Int_t lNbBinsPerVar[3] = {0};
     lNbBinsPerVar[0] = 100;
     lNbBinsPerVar[1] = 1000;
     lNbBinsPerVar[2] = 22;
     fCFContCascadePIDAsOmegaMinus = new AliCFContainer(Form("fCFContCascadePIDAsOmegaMinus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
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
     fListHistCascade->Add(fCFContCascadePIDAsOmegaMinus);
  }
  // PID container Omega+
  if(! fCFContCascadePIDAsOmegaPlus)  {
     const Int_t  lNbSteps      =  7;
     const Int_t  lNbVariables  =  3;
       //Array for the number of bins in each dimension :
     Int_t lNbBinsPerVar[3]= {0};
     lNbBinsPerVar[0] = 100;
     lNbBinsPerVar[1] = 1000;
     lNbBinsPerVar[2] = 22;  
     fCFContCascadePIDAsOmegaPlus = new AliCFContainer(Form("fCFContCascadePIDAsOmegaPlus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
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
     fListHistCascade->Add(fCFContCascadePIDAsOmegaPlus);
  }
  // Container for optimisation of topological selections 
  if(! fCFContAsCascadeCuts){
	// Container meant to store all the relevant distributions corresponding to the cut variables.
        //          - NB overflow/underflow of variables on which we want to cut later should be 0!!!
     const Int_t  lNbSteps      =  4;
     const Int_t  lNbVariables  =  19;
       //Array for the number of bins in each dimension :
     Int_t lNbBinsPerVar[lNbVariables] = {0};
     lNbBinsPerVar[0]  = 25;   //DcaCascDaughters                : [0.0,2.,3.0]        -> Rec.Cut = 2.0; 
     lNbBinsPerVar[1]  = 25;   //DcaBachToPrimVertex             : [0.0,0.24,100.0]    -> Rec.Cur = 0.01;
     lNbBinsPerVar[2]  = 30;   //CascCosineOfPointingAngle       : [0.97,1.]           -> Rec.Cut = 0.98;
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
     fCFContAsCascadeCuts = new AliCFContainer(Form("fCFContAsCascadeCuts_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Cut Container for Asso. Cascades", lNbSteps, lNbVariables, lNbBinsPerVar );
       //Setting the bin limits 
       //0 - DcaCascDaughters
     Double_t *lBinLim0 = new Double_t[ lNbBinsPerVar[0]+1 ];
     for(Int_t i=0; i<lNbBinsPerVar[0]; i++) lBinLim0[i] = (Double_t)0.0 + (2.4 -0.0)/(lNbBinsPerVar[0] - 1) * (Double_t)i;
     lBinLim0[ lNbBinsPerVar[0] ] = 3.0;
     fCFContAsCascadeCuts -> SetBinLimits(0, lBinLim0);
     delete[] lBinLim0;
       //1 - DcaBachToPrimVertex
     Double_t *lBinLim1 = new Double_t[ lNbBinsPerVar[1]+1 ];
     for(Int_t i=0; i<lNbBinsPerVar[1]; i++) lBinLim1[i] = (Double_t)0.0 + (0.24 - 0.0)/(lNbBinsPerVar[1] - 1) * (Double_t)i;
     lBinLim1[ lNbBinsPerVar[1] ] = 100.0;
     fCFContAsCascadeCuts -> SetBinLimits(1, lBinLim1);
     delete [] lBinLim1;
       //2 - CascCosineOfPointingAngle
     fCFContAsCascadeCuts -> SetBinLimits(2, .97, 1.);        
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
     fCFContAsCascadeCuts->SetVarTitle(6,  "cos(V0 PA) in cascade");
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
     fListHistCascade->Add(fCFContAsCascadeCuts);
  }

 PostData(1, fListHistCascade); 
 PostData(2, fCFContCascadePIDAsXiMinus);
 PostData(3, fCFContCascadePIDAsXiPlus);
 PostData(4, fCFContCascadePIDAsOmegaMinus);
 PostData(5, fCFContCascadePIDAsOmegaPlus);
 PostData(6, fCFContAsCascadeCuts);

}// end CreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskCheckPerformanceCascadepp276::UserExec(Option_t *) {
	
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


  //////////////////	
  // Event selection	
  //////////////////
  // In order:
  // 1) SDD selection
  // 2) Physics selection
  // 3) Select only looking at events with well-established PV
  // 4) Pileup selection
  // 5) |Z| < 10 cm
    
   //---------------------------------------------------------
   // Load the InputEvent and check it (for the ESD and AOD)
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
       // - Cascade vertexer (ESD)
       // Relaunch V0 and Cascade vertexer
       if (fkRerunV0CascVertexers) { 
           lESDevent->ResetCascades();
           lESDevent->ResetV0s();
           //AliV0vertexer lV0vtxer;
           //AliCascadeVertexer lCascVtxer;
           //lV0vtxer.SetCuts(fV0Sels);
           //lCascVtxer.SetCuts(fCascSels);
           //lV0vtxer.Tracks2V0vertices(lESDevent);
           //lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
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

   //------------------------------
   // - Plots Before any selections
   //------------------------------
   // - Define variables
   Int_t  ncascadesBeforeAnySel          = -1; //number of cascades before any selections
   Int_t  nTrackMultiplicityBeforeAnySel = -1; //number of tracks before any selections
   if (fAnalysisType == "ESD") {
       //Multiplicity
       Int_t lMultiplicity = -100;
       lMultiplicity = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5);
       nTrackMultiplicityBeforeAnySel = lMultiplicity;
       ncascadesBeforeAnySel = lESDevent->GetNumberOfCascades();
   } else if (fAnalysisType == "AOD") {
       //Multiplicity
       Int_t lMultiplicity = -100;
       nTrackMultiplicityBeforeAnySel = lMultiplicity;
       ncascadesBeforeAnySel  = lAODevent->GetNumberOfCascades();
   }
   fHistTrackMultiplicityBeforeAnySel->Fill(nTrackMultiplicityBeforeAnySel);
   fHistCascadeMultiplicityBeforeAnySel->Fill(ncascadesBeforeAnySel);

   //----------------
   // - SDD selection
   //----------------
   // - Define variables
   Int_t  ncascadesAfterSDDSel          = -1; //number of cascades after SDD selection
   Int_t  nTrackMultiplicityAfterSDDSel = -1; //number of tracks after SDD selection
   if (fkSDDselectionOn) {
        TString trcl = " ";
        trcl = lESDevent->GetFiredTriggerClasses();
        if (fAnalysisType == "ESD") trcl = lESDevent->GetFiredTriggerClasses();
        else if (fAnalysisType == "AOD") trcl = lAODevent->GetFiredTriggerClasses();
        if (fwithSDD){   // ---> Select event with SDD ON
            if(!(trcl.Contains("ALLNOTRD"))) {
                 PostData(1, fListHistCascade);
                 PostData(2, fCFContCascadePIDAsXiMinus);
                 PostData(3, fCFContCascadePIDAsXiPlus);
                 PostData(4, fCFContCascadePIDAsOmegaMinus);
                 PostData(5, fCFContCascadePIDAsOmegaPlus);
                 PostData(6, fCFContAsCascadeCuts);
                 cout<<"Bad event: SDD turn OFF =>  RETURN!! (Exclude it)..."<<endl;
                 return;
            } else {
                 cout<<"Good event: SDD turn ON."<<endl;
            }
        } else if (!fwithSDD){  // ---> Select event with SDD OFF
            if((trcl.Contains("ALLNOTRD"))) {
                 PostData(1, fListHistCascade);
                 PostData(2, fCFContCascadePIDAsXiMinus);
                 PostData(3, fCFContCascadePIDAsXiPlus);
                 PostData(4, fCFContCascadePIDAsOmegaMinus);
                 PostData(5, fCFContCascadePIDAsOmegaPlus);
                 PostData(6, fCFContAsCascadeCuts);
                 cout<<"Bad event:  SDD turn ON =>  RETURN!! (Exclude it)..."<<endl;
                 return;
            } else {
                 cout<<"Good event: SDD turn OFF."<<endl;
            }
        }
   }
   // - Take the number of cascades and tracks after the SDD selection
   if (fAnalysisType == "ESD") {
       Int_t lMultiplicity = -100;
       lMultiplicity = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5);
       ncascadesAfterSDDSel = lESDevent->GetNumberOfCascades();
       nTrackMultiplicityAfterSDDSel = lMultiplicity;
   } else if (fAnalysisType == "AOD") {
       Int_t lMultiplicity = -100;
       ncascadesAfterSDDSel = lAODevent->GetNumberOfCascades();
       nTrackMultiplicityAfterSDDSel = lMultiplicity;
   }
   // - Fill the plots
   fHistTrackMultiplicityAfterSDDSel->Fill(nTrackMultiplicityAfterSDDSel);
   fHistCascadeMultiplicityAfterSDDSel->Fill(ncascadesAfterSDDSel);

   //------------------------------
   // - Plots pre-physics selection
   //------------------------------
   // - Produce the 3Dhisto for the efficiency denominator
   Int_t lNbMCPrimary = 0;
   lNbMCPrimary = lMCstack->GetNprimary();

   for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < lNbMCPrimary; iCurrentLabelStack++) {

     Double_t partEnergy = 0.;
     Double_t partPz     = 0.;
     Double_t partP      = 0.;
     Double_t partPt     = 0.;
     Double_t partVx     = 0.;
     Double_t partVy     = 0.;
     Double_t partVz     = 0.;
     Double_t bacVx      = 0.;
     Double_t bacVy      = 0.;
     Double_t bacVz      = 0.;
     Double_t partMass   = 0.;
     Int_t    PDGcode    = 0;
     Int_t lPrimaryTrackMultiplicity = nTrackMultiplicityAfterSDDSel;

       if ( fAnalysisType == "ESD" ) {
            TParticle* lCurrentParticlePrimary = 0x0;
            lCurrentParticlePrimary = lMCstack->Particle( iCurrentLabelStack );        
            if (!lCurrentParticlePrimary) {
                  Printf("Cascade loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
                  continue;
            }
            if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) continue;
            TParticle* xiMC = 0x0;
            xiMC = lCurrentParticlePrimary;
            if (!xiMC) {
                  Printf("MC TParticle pointer to Cascade = 0x0 ! Skip ...");
                  continue;
            }
            partEnergy = xiMC->Energy();
            partPz     = xiMC->Pz();
            partPt     = xiMC->Pt();
            partP      = xiMC->P();
            partMass   = xiMC->GetMass();
            partVx     = xiMC->Vx();
            partVy     = xiMC->Vy();
            partVz     = xiMC->Vz();
            if (xiMC->GetDaughter(0)>=0) {    
                 TParticle *mcBach = lMCstack->Particle(xiMC->GetDaughter(0));
                 if (mcBach) {
                     bacVx  = mcBach->Vx();
                     bacVy  = mcBach->Vy();
                     bacVz  = mcBach->Vz();
                 }
            }
            PDGcode = lCurrentParticlePrimary->GetPdgCode();
       } else if ( fAnalysisType == "AOD" ) {
            AliAODMCParticle *lCurrentParticleaod = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack);
            if (!lCurrentParticleaod) {
                  Printf("Cascade loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
                  continue;
            }
            if (!lCurrentParticleaod->IsPhysicalPrimary()) continue;
            partEnergy = lCurrentParticleaod->E();
            partPz     = lCurrentParticleaod->Pz();
            partP      = lCurrentParticleaod->P();
            partPt     = lCurrentParticleaod->Pt();
            partMass   = lCurrentParticleaod->M(); 
            partVx     = lCurrentParticleaod->Xv();
            partVy     = lCurrentParticleaod->Yv();
            partVz     = lCurrentParticleaod->Zv();
            if (lCurrentParticleaod->GetDaughter(0)>=0) {
                 AliAODMCParticle *mcBach = (AliAODMCParticle*) arrayMC->At(lCurrentParticleaod->GetDaughter(0));
                 if (mcBach) {
                     bacVx  = mcBach->Xv();
                     bacVy  = mcBach->Yv();
                     bacVz  = mcBach->Zv();
                 }
            }     
            PDGcode = lCurrentParticleaod->GetPdgCode();
       }

       // - Calculate rapidity
       Double_t lRapXiMC = 0.5*TMath::Log((partEnergy + partPz) / (partEnergy - partPz + 1.e-13));
       // - Calculate proper lenght
       Double_t lctau = TMath::Sqrt((partVx-bacVx)*(partVx-bacVx)+(partVy-bacVy)*(partVy-bacVy)+(partVz-bacVz)*(partVz-bacVz));
       if (partP != 0.) lctau = lctau*partMass/partP;
       else lctau = -1.;
       // - Fill Histograms
       if (PDGcode ==  3312) {
           f3dHistGenPtVsGenYvsNtracksXiMinus->Fill(partPt, lRapXiMC, lPrimaryTrackMultiplicity);
           f3dHistGenPtVsGenctauvsYXiMinus->Fill(partPt, lctau, lRapXiMC);
       }
       if (PDGcode == -3312) {
           f3dHistGenPtVsGenYvsNtracksXiPlus->Fill(partPt, lRapXiMC, lPrimaryTrackMultiplicity);
           f3dHistGenPtVsGenctauvsYXiPlus->Fill(partPt, lctau, lRapXiMC);
       }
       if (PDGcode ==  3334) {
           f3dHistGenPtVsGenYvsNtracksOmegaMinus->Fill(partPt, lRapXiMC, lPrimaryTrackMultiplicity);
           f3dHistGenPtVsGenctauvsYOmegaMinus->Fill(partPt, lctau, lRapXiMC);
       }
       if (PDGcode == -3334) {
           f3dHistGenPtVsGenYvsNtracksOmegaPlus->Fill(partPt, lRapXiMC, lPrimaryTrackMultiplicity);
           f3dHistGenPtVsGenctauvsYOmegaPlus->Fill(partPt, lctau, lRapXiMC);
       }
   }

 
   //--------------------
   // - Physics selection
   //--------------------
   // - Define new variables
   Int_t    ncascadesAfterPhysicsSel          = -1; //number of cascades after physics selection
   Int_t    nTrackMultiplicityAfterPhysicsSel = -1; //number of tracks after physics selection
   // - Selection for ESD and AOD
   if (fAnalysisType == "ESD") {
       UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
       Bool_t isSelected = 0;
       isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
       if(!isSelected){
           PostData(1, fListHistCascade);
           PostData(2, fCFContCascadePIDAsXiMinus);
           PostData(3, fCFContCascadePIDAsXiPlus);
           PostData(4, fCFContCascadePIDAsOmegaMinus);
           PostData(5, fCFContCascadePIDAsOmegaPlus);
           PostData(6, fCFContAsCascadeCuts);
           return;
       }
       // - Take the number of cascades and tracks after physics selection
       ncascadesAfterPhysicsSel = lESDevent->GetNumberOfCascades();
       nTrackMultiplicityAfterPhysicsSel = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5);
   } else if (fAnalysisType == "AOD") {
       UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
       Bool_t isSelected = 0;
       isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
       if(!isSelected){
           PostData(1, fListHistCascade);
           PostData(2, fCFContCascadePIDAsXiMinus);
           PostData(3, fCFContCascadePIDAsXiPlus);
           PostData(4, fCFContCascadePIDAsOmegaMinus);
           PostData(5, fCFContCascadePIDAsOmegaPlus);
           PostData(6, fCFContAsCascadeCuts);
           return;
       }
       // - Take the number of cascades and tracks after physics selection
       ncascadesAfterPhysicsSel = lAODevent->GetNumberOfCascades();
       nTrackMultiplicityAfterPhysicsSel = -100;
   }
   fHistCascadeMultiplicityAfterPhysicsSel->Fill(ncascadesAfterPhysicsSel);
   fHistTrackMultiplicityAfterPhysicsSel->Fill(nTrackMultiplicityAfterPhysicsSel);

   //-------------------------------------------------------
   // Select only looking at events with well-established PV
   //-------------------------------------------------------
   Int_t    ncascadesForSelEvtNoTPCOnly                  = -1; //number of cascades after the TPConly selection
   Int_t    nTrackMultiplicityForSelEvtNoTPCOnly         = -1; //number of tracks after the TPConly selection
   if (fAnalysisType == "ESD" ) {
       // - Select only looking at events with well-established PV
       if (fkQualityCutNoTPConlyPrimVtx) {
           const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
           const AliESDVertex *lPrimaryTrackingVtx = lESDevent->GetPrimaryVertexTracks();
           if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingVtx->GetStatus() ){
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
       // - Take the number of cascades and tracks after TPConly selection
       ncascadesForSelEvtNoTPCOnly = lESDevent->GetNumberOfCascades();
       nTrackMultiplicityForSelEvtNoTPCOnly = fESDtrackCuts->GetReferenceMultiplicity(lESDevent,AliESDtrackCuts::kTrackletsITSTPC,0.5);
   } else if (fAnalysisType == "AOD") {
       // - Select only looking at events with well-established PV
       if (fkQualityCutNoTPConlyPrimVtx) {
           const AliAODVertex *lPrimarySPDVtx = lAODevent->GetPrimaryVertexSPD();
           const AliAODVertex *lPrimaryTrackingAODVtx = lAODevent->GetPrimaryVertex();
           if (!lPrimarySPDVtx && !lPrimaryTrackingAODVtx) {
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
       // - Take the number of cascades and tracks after TPConly selection
       ncascadesForSelEvtNoTPCOnly = lAODevent->GetNumberOfCascades();
       nTrackMultiplicityForSelEvtNoTPCOnly = -100;  //FIXME
   }
   fHistCascadeMultiplicityForSelEvtNoTPCOnly->Fill(ncascadesForSelEvtNoTPCOnly);
   fHistTrackMultiplicityForSelEvtNoTPCOnly->Fill(nTrackMultiplicityForSelEvtNoTPCOnly);
    
   //-----------------
   // Pileup selection
   //-----------------
   Int_t    ncascadesForSelEvtNoTPCOnlyNoPileup          = -1; //number of cascades after the NoPileup selection
   Int_t    nTrackMultiplicityForSelEvtNoTPCOnlyNoPileup = -1; //number of tracks after the Pileup selection
   if (fAnalysisType == "ESD" ) {
       // - Selection for pile up
       if (fkRejectEventPileUp) {
           if(lESDevent->IsPileupFromSPD()){
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
       // - Take the number of cascades and tracks after Pileup selection
       ncascadesForSelEvtNoTPCOnlyNoPileup = lESDevent->GetNumberOfCascades();
       nTrackMultiplicityForSelEvtNoTPCOnlyNoPileup = fESDtrackCuts->GetReferenceMultiplicity(lESDevent,AliESDtrackCuts::kTrackletsITSTPC,0.5);
   } else if (fAnalysisType == "AOD") {
       // - Selection for pile up
       if (fkRejectEventPileUp) {
           if(lAODevent->IsPileupFromSPD()){
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
       // - Take the number of cascades and tracks after Pileup selection
       ncascadesForSelEvtNoTPCOnlyNoPileup = lAODevent->GetNumberOfCascades();
       nTrackMultiplicityForSelEvtNoTPCOnlyNoPileup = -100;
   }
   fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup->Fill(ncascadesForSelEvtNoTPCOnlyNoPileup);
   fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup->Fill(nTrackMultiplicityForSelEvtNoTPCOnlyNoPileup);
    
   //-------------------
   // - Vertex selection 
   //-------------------
   Int_t    ncascadesAfterVertexSel                      = -1; //number of cascades after vertex selection
   Int_t    nTrackMultiplicityAfterVertexSel             = -1; //number of tracks after vertex selection
   Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
   Double_t tPrimaryVtxPosition[3] = {-100.0, -100.0, -100.0}; 
   Double_t lMagneticField = -10.;
   if (fAnalysisType == "ESD" ) {
       // - Primary vertex definition
       const AliESDVertex *lPrimaryBestVtx = lESDevent->GetPrimaryVertex();
       if (!lPrimaryBestVtx) {
            AliWarning("No prim. vertex in AOD... return!");
            PostData(1, fListHistCascade);
            PostData(2, fCFContCascadePIDAsXiMinus);
            PostData(3, fCFContCascadePIDAsXiPlus);
            PostData(4, fCFContCascadePIDAsOmegaMinus);
            PostData(5, fCFContCascadePIDAsOmegaPlus);
            PostData(6, fCFContAsCascadeCuts);
            return;
       }
       lPrimaryBestVtx->GetXYZ( lBestPrimaryVtxPos );
       // - Vertex position before any event selection on vertex position
       const AliVVertex *primaryVtx = lESDevent->GetPrimaryVertex();
       tPrimaryVtxPosition[0] = primaryVtx->GetX();
       tPrimaryVtxPosition[1] = primaryVtx->GetY();
       tPrimaryVtxPosition[2] = primaryVtx->GetZ();
       fHistPVx->Fill( tPrimaryVtxPosition[0] );
       fHistPVy->Fill( tPrimaryVtxPosition[1] );
       fHistPVz->Fill( tPrimaryVtxPosition[2] );
       // - Get magnetic filed info
       lMagneticField = lESDevent->GetMagneticField();
       // - Selection on the primary vertex Z position 
       if (fkQualityCutZprimVtxPos) {
           if (TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRange || TMath::Abs(lBestPrimaryVtxPos[2]) < fVtxRangeMin) {
                AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !");
                PostData(1, fListHistCascade);
                PostData(2, fCFContCascadePIDAsXiMinus);
                PostData(3, fCFContCascadePIDAsXiPlus);
                PostData(4, fCFContCascadePIDAsOmegaMinus);
                PostData(5, fCFContCascadePIDAsOmegaPlus);
                PostData(6, fCFContAsCascadeCuts);
                return;
          }
       }
       // - Take the number of cascades and tracks after vertex Z position selection
       ncascadesAfterVertexSel = lESDevent->GetNumberOfCascades();
       nTrackMultiplicityAfterVertexSel = fESDtrackCuts->GetReferenceMultiplicity(lESDevent,AliESDtrackCuts::kTrackletsITSTPC,0.5);
   } else if (fAnalysisType == "AOD") {
       // - Primary vertex definition
       const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex(); // get the best primary vertex available for the event GetVertex(0)
       if (!lPrimaryBestAODVtx) {
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
       // - Vertex position before any event selection on vertex position
       const AliVVertex *primaryVtx = lAODevent->GetPrimaryVertex();
       tPrimaryVtxPosition[0] = primaryVtx->GetX();
       tPrimaryVtxPosition[1] = primaryVtx->GetY();
       tPrimaryVtxPosition[2] = primaryVtx->GetZ();
       fHistPVx->Fill( tPrimaryVtxPosition[0] );
       fHistPVy->Fill( tPrimaryVtxPosition[1] );
       fHistPVz->Fill( tPrimaryVtxPosition[2] );
       // - Get magnetic filed info
       lMagneticField = lAODevent->GetMagneticField();
       // - Selection on the primary vertex Z position 
       if (fkQualityCutZprimVtxPos) {
           if (TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRange && TMath::Abs(lBestPrimaryVtxPos[2]) < fVtxRangeMin) {
                AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !");
                PostData(1, fListHistCascade);
                PostData(2, fCFContCascadePIDAsXiMinus);
                PostData(3, fCFContCascadePIDAsXiPlus);
                PostData(4, fCFContCascadePIDAsOmegaMinus);
                PostData(5, fCFContCascadePIDAsOmegaPlus);
                PostData(6, fCFContAsCascadeCuts);
                return;
           }
       }
       // - Take the number of cascades and tracks after vertex Z position selection
       ncascadesAfterVertexSel = lAODevent->GetNumberOfCascades();
       nTrackMultiplicityAfterVertexSel = -100;
   }
   // - Fill the plots
   fHistCascadeMultiplicityAfterVertexCutSel->Fill(ncascadesAfterVertexSel);
   fHistTrackMultiplicityAfterVertexCutSel->Fill(nTrackMultiplicityAfterVertexSel);

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
    

   //----------------------------------------------------------------------	
   // - Loop over the different types of GENERATED cascades (Xi-+, Omega-+)	
   //----------------------------------------------------------------------
   // - Initialisation of useful local variables		
   Int_t lPdgCodeCasc            = 0;
   Int_t lPdgCodeBach            = 0;
   Int_t lPdgCodeLambda          = 0;
   Int_t lPdgCodeDghtMesV0       = 0;
   Int_t lPdgCodeDghtBarV0       = 0;	
   TH1F *lHistEtaGenCasc         = 0;	
   TH3D *l3dHistGenPtVsGenYvsNtracksPhysEff = 0;
   TH3D *l3dHistGenPtVsGenctauvsYPhysEff    = 0;
   TH1F *lHistThetaGenCasc       = 0;
   TH2D *l2dHistGenPtVsGenYFdbl  = 0;
   TH1F *lHistThetaLambda        = 0;
   TH1F *lHistThetaBach          = 0;
   TH1F *lHistThetaBarDghter     = 0;
   TH1F *lHistThetaMesDghter     = 0;
   TH1F *lHistPtBach             = 0;
   TH1F *lHistPtBarDghter        = 0;
   TH1F *lHistPtMesDghter        = 0;
   Int_t ncascperev = 0; 
   Int_t ncascperevtot = 0;

   for (Int_t iCascType = 1; iCascType < 5; iCascType++) { 
         ncascperev = 0;
         ncascperevtot = 0;
         Int_t lPrimaryTrackMultiplicity = nTrackMultiplicityAfterSDDSel;

         switch (iCascType) {
           case 1: // Xi-
               lPdgCodeCasc       =   3312;  //Xi-
               lPdgCodeBach       =   -211;  //Pi-
               lPdgCodeLambda     =   3122;  //Lambda0
               lPdgCodeDghtMesV0  =   -211;  //Pi-
               lPdgCodeDghtBarV0  =   2212;  //Proton 	
               lHistEtaGenCasc        = fHistEtaGenCascXiMinus;         // this plot for any Xi- 
  	       lHistThetaGenCasc      = fHistThetaGenCascXiMinus;       // cascades generated within acceptance (cut in pt + theta)
               l3dHistGenPtVsGenYvsNtracksPhysEff = f3dHistGenPtVsGenYvsNtracksXiMinusPhysEff;
               l3dHistGenPtVsGenctauvsYPhysEff    = f3dHistGenPtVsGenctauvsYXiMinusPhysEff;
	       l2dHistGenPtVsGenYFdbl = f2dHistGenPtVsGenYFdblXiMinus;
	       lHistThetaLambda       = fHistThetaLambdaXiMinus;
	       lHistThetaBach         = fHistThetaBachXiMinus;
	       lHistThetaBarDghter    = fHistThetaBarDghterXiMinus;
	       lHistThetaMesDghter    = fHistThetaMesDghterXiMinus;
	       lHistPtBach	      = fHistPtBachXiMinus;
	       lHistPtBarDghter       = fHistPtBarDghterXiMinus;
	       lHistPtMesDghter       = fHistPtMesDghterXiMinus;
               break; 
           case 2: // Xi+
               lPdgCodeCasc        =  -3312;  //Xi+
               lPdgCodeBach        =    211;  //Pi+
               lPdgCodeLambda      =  -3122;  //AntiLambda0
               lPdgCodeDghtMesV0   =    211;  //Pi+
               lPdgCodeDghtBarV0   =  -2212;  //AntiProton  
      	       lHistEtaGenCasc        = fHistEtaGenCascXiPlus;       // this plot for any Xi+
	       lHistThetaGenCasc      = fHistThetaGenCascXiPlus;     // cascades generated within acceptance (cut in pt + theta)
               l3dHistGenPtVsGenYvsNtracksPhysEff = f3dHistGenPtVsGenYvsNtracksXiPlusPhysEff;
               l3dHistGenPtVsGenctauvsYPhysEff    = f3dHistGenPtVsGenctauvsYXiPlusPhysEff;
	       l2dHistGenPtVsGenYFdbl = f2dHistGenPtVsGenYFdblXiPlus;
	       lHistThetaLambda       = fHistThetaLambdaXiPlus;
	       lHistThetaBach         = fHistThetaBachXiPlus;
	       lHistThetaBarDghter    = fHistThetaBarDghterXiPlus;
	       lHistThetaMesDghter    = fHistThetaMesDghterXiPlus;
	       lHistPtBach	      = fHistPtBachXiPlus;
	       lHistPtBarDghter       = fHistPtBarDghterXiPlus;
	       lHistPtMesDghter       = fHistPtMesDghterXiPlus;  
    	       break;
           case 3: // Omega-
    	       lPdgCodeCasc       =   3334;  //Omega-
               lPdgCodeBach       =   -321;  //K-
               lPdgCodeLambda     =   3122;  //Lambda0
               lPdgCodeDghtMesV0  =   -211;  //Pi-
               lPdgCodeDghtBarV0  =   2212;  //Proton
	       lHistEtaGenCasc        = fHistEtaGenCascOmegaMinus;        // this plot for any Omega+	 	
	       lHistThetaGenCasc      = fHistThetaGenCascOmegaMinus;      // cascades generated within acceptance (cut in pt + theta)
	       l2dHistGenPtVsGenYFdbl = f2dHistGenPtVsGenYFdblOmegaMinus;
               l3dHistGenPtVsGenYvsNtracksPhysEff = f3dHistGenPtVsGenYvsNtracksOmegaMinusPhysEff;
               l3dHistGenPtVsGenctauvsYPhysEff    = f3dHistGenPtVsGenctauvsYOmegaMinusPhysEff;
	       lHistThetaLambda       = fHistThetaLambdaOmegaMinus;
	       lHistThetaBach         = fHistThetaBachOmegaMinus;
	       lHistThetaBarDghter    = fHistThetaBarDghterOmegaMinus;
	       lHistThetaMesDghter    = fHistThetaMesDghterOmegaMinus;
	       lHistPtBach	      = fHistPtBachOmegaMinus;
	       lHistPtBarDghter       = fHistPtBarDghterOmegaMinus;
	       lHistPtMesDghter       = fHistPtMesDghterOmegaMinus;   
               break;
           case 4:  // Omega+
               lPdgCodeCasc       =  -3334;  //Omega+
               lPdgCodeBach       =    321;  //K+
               lPdgCodeLambda     =  -3122;  //AntiLambda0
               lPdgCodeDghtMesV0  =    211;  //Pi+
               lPdgCodeDghtBarV0  =  -2212;  //AntiProton 
	       lHistEtaGenCasc        = fHistEtaGenCascOmegaPlus;        // this plot for any Omega-
	       lHistThetaGenCasc      = fHistThetaGenCascOmegaPlus;      // cascades generated within acceptance (cut in pt + theta)
	       l2dHistGenPtVsGenYFdbl = f2dHistGenPtVsGenYFdblOmegaPlus;
               l3dHistGenPtVsGenYvsNtracksPhysEff = f3dHistGenPtVsGenYvsNtracksOmegaPlusPhysEff;
               l3dHistGenPtVsGenctauvsYPhysEff    = f3dHistGenPtVsGenctauvsYOmegaPlusPhysEff;
	       lHistThetaLambda       = fHistThetaLambdaOmegaPlus;
	       lHistThetaBach         = fHistThetaBachOmegaPlus;
	       lHistThetaBarDghter    = fHistThetaBarDghterOmegaPlus;
	       lHistThetaMesDghter    = fHistThetaMesDghterOmegaPlus;
	       lHistPtBach	      = fHistPtBachOmegaPlus;
	       lHistPtBarDghter       = fHistPtBarDghterOmegaPlus;
	       lHistPtMesDghter       = fHistPtMesDghterOmegaPlus;  
               break;
         }

         for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < lNbMCPrimary; iCurrentLabelStack++) {

               Double_t partEnergy = 0.;
               Double_t partPz     = 0.;
               Double_t partEta    = 0.;
               Double_t partTheta  = 0.;
               Double_t partP      = 0.;
               Double_t partPt     = 0.;
               Double_t partVx     = 0.;
               Double_t partVy     = 0.; 
               Double_t partVz     = 0.;
               Double_t bacVx      = 0.;
               Double_t bacVy      = 0.;
               Double_t bacVz      = 0.;    
               Double_t partMass   = 0.;

               if ( fAnalysisType == "ESD" ) {      
                    TParticle* lCurrentParticle = 0x0; 
                    lCurrentParticle = lMCstack->Particle( iCurrentLabelStack );
                    if (!lCurrentParticle) {
                        Printf("Cascade loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
                        continue;
                    }
                    if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) continue; 
                    if (lCurrentParticle->GetPdgCode() == lPdgCodeCasc) {  // Here !	
	                TParticle* xiMC = 0x0;
	                xiMC = lCurrentParticle;
	                if (!xiMC) {
	                    Printf("MC TParticle pointer to Cascade = 0x0 ! Skip ...");
	                    continue;
	                }
                        partEnergy = xiMC->Energy();
                        partPz     = xiMC->Pz();
                        partEta    = xiMC->Eta();
                        partPt     = xiMC->Pt();
                        partP      = xiMC->P();
                        partTheta  = xiMC->Theta();
                        partMass   = xiMC->GetMass();
                        partVx     = xiMC->Vx();
                        partVy     = xiMC->Vy();
                        partVz     = xiMC->Vz();
                        if (xiMC->GetDaughter(0)>=0) {
                             TParticle *mcBach = lMCstack->Particle(xiMC->GetDaughter(0));
                             if (mcBach) {
                                 bacVx  = mcBach->Vx();
                                 bacVy  = mcBach->Vy();
                                 bacVz  = mcBach->Vz();
                             }
                        }
                    } else continue;
               } else if ( fAnalysisType == "AOD" ) {
                    AliAODMCParticle *lCurrentParticleaod = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack);
                    if (!lCurrentParticleaod) {
                        Printf("Cascade loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
                        continue;
                    }
                    if (!lCurrentParticleaod->IsPhysicalPrimary()) continue;  
                    if (!(lCurrentParticleaod->PdgCode() == lPdgCodeCasc)) continue;
                    partEnergy = lCurrentParticleaod->E();
                    partPz     = lCurrentParticleaod->Pz();
                    partEta    = lCurrentParticleaod->Eta();
                    partP      = lCurrentParticleaod->P();
                    partPt     = lCurrentParticleaod->Pt();
                    partTheta  = lCurrentParticleaod->Theta();
                    partMass   = lCurrentParticleaod->M();   //FIXME: not sure this works, seems not implemented
                    partVx     = lCurrentParticleaod->Xv();
                    partVy     = lCurrentParticleaod->Yv();
                    partVz     = lCurrentParticleaod->Zv();
                    if (lCurrentParticleaod->GetDaughter(0)>=0) {
                         AliAODMCParticle *mcBach = (AliAODMCParticle*) arrayMC->At(lCurrentParticleaod->GetDaughter(0));
                         if (mcBach) {
                              bacVx  = mcBach->Xv();
                              bacVy  = mcBach->Yv();
                              bacVz  = mcBach->Zv();
                         } 
                    }
               }
               ncascperevtot++;	
               // - Fill the first histos : = any generated Xi, not necessarily within the acceptance
               Double_t lRapXiMC = 0.5*TMath::Log((partEnergy + partPz) / (partEnergy - partPz +1.e-13));
               // - Calculate proper time
               Double_t lctau = TMath::Sqrt((partVx-bacVx)*(partVx-bacVx)+(partVy-bacVy)*(partVy-bacVy)+(partVz-bacVz)*(partVz-bacVz));
               if (partP!=0.)    lctau = lctau*partMass/partP;
               else lctau = -1.;
               Double_t lRadToDeg = 180.0/TMath::Pi();
               // - Fill the first histos : = any generated Xi, not necessarily within the acceptance		
               lHistEtaGenCasc->Fill( partEta );	 
               l3dHistGenPtVsGenYvsNtracksPhysEff->Fill( partPt, lRapXiMC, lPrimaryTrackMultiplicity );
               l3dHistGenPtVsGenctauvsYPhysEff->Fill( partPt, lctau, lRapXiMC );
               lHistThetaGenCasc->Fill( lRadToDeg * partTheta );

               //--------------------------------------------------------------------------------------------
               // - Check the emission of particle stays within the acceptance of the detector (cut in theta)
               if (fApplyAccCut) { if( partTheta < TMath::Pi()/4.0 || partTheta > 3.0*TMath::Pi()/4.0 ) continue;}	

               Float_t lambdaTheta = 0.;
               Float_t bacTheta    = 0.;
               Float_t dghtBarV0Theta = 0.;
               Float_t dghtMesV0Theta = 0.;
               Float_t bacPt       = 0.;
               Float_t dghtBarV0Pt = 0.;
               Float_t dghtMesV0Pt = 0.;

               if ( fAnalysisType == "ESD" ) { 
                    TParticle* xiMC = lMCstack->Particle( iCurrentLabelStack ); 
                    if ( xiMC->GetNDaughters() != 2) continue;
                    if ( xiMC->GetDaughter(0) < 0 )  continue;
                    if ( xiMC->GetDaughter(1) < 0 )  continue;	
                    TParticle* lDght0ofXi = lMCstack->Particle(  xiMC->GetDaughter(0) );
                    TParticle* lDght1ofXi = lMCstack->Particle(  xiMC->GetDaughter(1) );
                    TParticle* lLambda = 0;
                    TParticle* lBach   = 0;

                    // Xi - Case 1
                    if ( lDght0ofXi->GetPdgCode() == lPdgCodeLambda && lDght1ofXi->GetPdgCode() == lPdgCodeBach ){  			
                           lLambda = lDght0ofXi;   // dghter0 = Lambda
                           lBach   = lDght1ofXi;   // dghter1 = Pi-
                    }  		
                    // Xi - Case 2
                    else if ( lDght0ofXi->GetPdgCode() == lPdgCodeBach && lDght1ofXi->GetPdgCode() == lPdgCodeLambda ){ 
                           lBach   = lDght0ofXi;   // dghter0 = Pi-
	                   lLambda = lDght1ofXi;   // dghter1 = Lambda
                    }
	            // Otherwise - Case 3	
                    else continue;

                    // - Check the emission of particle stays within the acceptance of the detector (cut in pt + theta)
                    if (fApplyAccCut) { 
                         if( lLambda->Theta() < TMath::Pi()/4.0 || lLambda->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
                         if( lBach->Theta() < TMath::Pi()/4.0   || lBach->Theta() > 3.0*TMath::Pi()/4.0 )   continue;
                         if( lBach->Pt() < 0.150 ) continue; //FIXME: maybe tuned for Xi but not for K- from Omega ...
                    } 				
 
                    //---------
                    // - V0 level
                    TParticle* lDghtBarV0 = 0;
                    TParticle* lDghtMesV0 = 0;
                    if( lLambda->GetNDaughters() != 2 )  continue;
                    if( lLambda->GetDaughter(0) < 0 )    continue;
                    if( lLambda->GetDaughter(1) < 0 )    continue;
                    TParticle* lDght0ofLambda = lMCstack->Particle(  lLambda->GetDaughter(0) );
                    TParticle* lDght1ofLambda = lMCstack->Particle(  lLambda->GetDaughter(1) );

                    // V0 - Case 1
                    if ( lDght0ofLambda->GetPdgCode() == lPdgCodeDghtBarV0 && lDght1ofLambda->GetPdgCode() == lPdgCodeDghtMesV0 ) {    // Here !			
                           lDghtBarV0 = lDght0ofLambda;   // dghter0 = Proton
                           lDghtMesV0 = lDght1ofLambda;   // dghter1 = Pi-
                    }  		
                    // V0 - Case 2
                    else if ( lDght0ofLambda->GetPdgCode() == lPdgCodeDghtMesV0 && lDght1ofLambda->GetPdgCode() == lPdgCodeDghtBarV0 ) {      // Here !
                           lDghtMesV0 = lDght0ofLambda;   // dghter0 = Pi-
	                   lDghtBarV0 = lDght1ofLambda;   // dghter1 = Proton
                    }		
                    // Otherwise - Case 3
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
                    bacPt          = lBach->Pt();
                    dghtBarV0Pt    = lDghtBarV0->Pt();
                    dghtMesV0Pt    = lDghtMesV0->Pt();
	  	 	
               } else if ( fAnalysisType == "AOD") {

                    AliAODMCParticle *xiMC = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack);
                    if (xiMC->GetNDaughters() != 2) continue;
                    if (xiMC->GetDaughter(0) < 0 )  continue;
                    if (xiMC->GetDaughter(1) < 0 )  continue;

                    AliAODMCParticle* lDght0ofXi = (AliAODMCParticle*) arrayMC->At(  xiMC->GetDaughter(0) );
                    AliAODMCParticle* lDght1ofXi = (AliAODMCParticle*) arrayMC->At(  xiMC->GetDaughter(1) );

                    AliAODMCParticle* lLambda = 0;
                    AliAODMCParticle* lBach   = 0;

                    // Xi - Case 1
                    if ( lDght0ofXi->PdgCode() == lPdgCodeLambda  &&  lDght1ofXi->PdgCode() == lPdgCodeBach ){ 
                            lLambda = lDght0ofXi;   // dghter0 = Lambda
                            lBach   = lDght1ofXi;   // dghter1 = Pi-
                    }
                    // Xi - Case 2
                    else if ( lDght0ofXi->PdgCode() == lPdgCodeBach && lDght1ofXi->PdgCode() == lPdgCodeLambda ){
                            lBach   = lDght0ofXi;   // dghter0 = Pi
                            lLambda = lDght1ofXi;   //dghter1 = Lambda
                    }
                    // Otherwise - Case 3
                    else continue;

                    // - Check the emission of particle stays within the acceptance of the detector (cut in pt + theta)
                    if (fApplyAccCut) {
                          if ( lLambda->Theta() < TMath::Pi()/4.0  ||    lLambda->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
                          if( lBach->Theta() < TMath::Pi()/4.0    ||    lBach->Theta() > 3.0*TMath::Pi()/4.0 )   continue;
                          if( lBach->Pt() < 0.150 ) continue; //FIXME : maybe tuned for Xi but not for K- from Omega ...
                    }

                    //-----------
                    // - V0 level 
                    AliAODMCParticle* lDghtBarV0 = 0;
                    AliAODMCParticle* lDghtMesV0 = 0;

                    if( lLambda->GetNDaughters() != 2 )  continue;
                    if( lLambda->GetDaughter(0) < 0 )    continue;
                    if( lLambda->GetDaughter(1) < 0 )    continue;

                    AliAODMCParticle* lDght0ofLambda = (AliAODMCParticle*) arrayMC->At(  lLambda->GetDaughter(0) );
                    AliAODMCParticle* lDght1ofLambda = (AliAODMCParticle*) arrayMC->At(  lLambda->GetDaughter(1) );

                    // V0 - Case 1
                    if ( lDght0ofLambda->PdgCode() == lPdgCodeDghtBarV0 && lDght1ofLambda->PdgCode() == lPdgCodeDghtMesV0 ) { 
                            lDghtBarV0 = lDght0ofLambda;   // dghter0 = Proton
                            lDghtMesV0 = lDght1ofLambda;   // dghter1 = Pi-
                    } 
                    // V0 - Case 2
                    else if ( lDght0ofLambda->PdgCode() == lPdgCodeDghtMesV0 && lDght1ofLambda->PdgCode() == lPdgCodeDghtBarV0 ) { 
                            lDghtMesV0 = lDght0ofLambda;   // dghter0 = Pi-
                            lDghtBarV0 = lDght1ofLambda;   // dghter1 = proton
                    } 
                    // V0 otherwise - Case 3
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
                    bacPt          = lBach->Pt();
                    dghtBarV0Pt    = lDghtBarV0->Pt();
                    dghtMesV0Pt    = lDghtMesV0->Pt();
               }

               //---------------------------------------
               // - Filling histos for findable cascades
               // - Fill theta histos 
               lHistThetaLambda->Fill( lRadToDeg * lambdaTheta );
               lHistThetaBach->Fill( lRadToDeg * bacTheta );
               lHistThetaBarDghter->Fill( lRadToDeg * dghtBarV0Theta );
               lHistThetaMesDghter->Fill( lRadToDeg * dghtMesV0Theta );
               // - Fill pt histos
               lHistPtBach             ->Fill( bacPt );
               lHistPtBarDghter        ->Fill( dghtBarV0Pt );
               lHistPtMesDghter        ->Fill( dghtMesV0Pt );
               l2dHistGenPtVsGenYFdbl  ->Fill( partPt, lRapXiMC );

               ncascperev++;			
			     
         }// This is the end of the loop on primaries
  
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
         lHistEtaGenCasc         = 0x0;
         lHistThetaGenCasc       = 0x0;
         l2dHistGenPtVsGenYFdbl  = 0x0;
         lHistThetaLambda        = 0x0;
         lHistThetaBach          = 0x0;
         lHistThetaBarDghter     = 0x0;
         lHistThetaMesDghter     = 0x0;
         lHistPtBach             = 0x0;
         lHistPtBarDghter        = 0x0;
         lHistPtMesDghter        = 0x0;	

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
   Double_t lDcaXiDaughters            = -1. ;
   Double_t lDcaBachToPrimVertexXi     = -1. ;
   Double_t lXiCosineOfPointingAngle   = -1. ;
   Double_t lPosXi[3]                  = { -1000.0, -1000.0, -1000.0 };
   Double_t lXiRadius                  = -1000. ;
   Double_t lInvMassLambdaAsCascDghter = 0.;
   Double_t lDcaV0DaughtersXi          = -1.;
   Double_t lV0CosineOfPointingAngleXi = -1.;
   Double_t lV0CosineOfPointingAngle   = -1.;
   Double_t lPosV0Xi[3]                = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
   Double_t lV0RadiusXi                = -1000.;
   Double_t lDcaV0ToPrimVertexXi       = -1.;
   Double_t lDcaPosToPrimVertexXi      = -1.;
   Double_t lDcaNegToPrimVertexXi      = -1.;
   Double_t lChargeXi                  = -1.;
   Double_t lV0mom                     = -1000.;
   Double_t lmcPt                      = -1.;         
   Double_t lmcRapCasc                 = -1.; 
   Double_t lmcEta                     = -1000.;     
   Double_t lmcTransvRadius            = -1000.; 
   Double_t lrecoPt                    = -100.;     
   Double_t lrecoTransvRadius          = -1000.; 
   Double_t lDeltaPhiMcReco            = -1.;
   Double_t lBachTransvMom             = 0.;
   Double_t lpTrackTransvMom           = 0.;
   Double_t lnTrackTransvMom           = 0.;
   Double_t lmcPtPosV0Dghter           = -100.;
   Double_t lmcPtNegV0Dghter           = -100.;
   Double_t lrecoP                     = -100.;
   Double_t lmcPtBach                  = -100.;
   Double_t cascadeMass                = 0.;

   // - Get the number of cascades
   Int_t ncascades = 0;
   if      ( fAnalysisType == "ESD" ) { ncascades = lESDevent->GetNumberOfCascades(); }
   else if ( fAnalysisType == "AOD" ) { ncascades = lAODevent->GetNumberOfCascades(); }

   //-------------------------------
   // - Begining of the Cascade Loop
   for (Int_t iXi = 0; iXi < ncascades; iXi++) {

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
        Bool_t   lIsBachelorKaonForTPC = kFALSE;
        Bool_t   lIsBachelorPionForTPC = kFALSE;
        Bool_t   lIsNegPionForTPC      = kFALSE;
        Bool_t   lIsPosPionForTPC      = kFALSE;
        Bool_t   lIsNegProtonForTPC    = kFALSE;
        Bool_t   lIsPosProtonForTPC    = kFALSE;

        // - Combined PID
        // Reasonable guess for the priors for the cascade track sample (e-, mu, pi, K, p)
        Double_t lPriorsGuessXi[14]    = {0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        Double_t lPriorsGuessOmega[14] = {0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        Double_t ppionBach = 0.0, pkaonBach = 0.0;
        Bool_t   lIsBachelorMCPiMinus  = kFALSE;
        Bool_t   lIsBachelorMCPiPlus   = kFALSE;
        Bool_t   lIsBachelorMCKMinus   = kFALSE;
        Bool_t   lIsBachelorMCKPlus    = kFALSE;
        Double_t lInvMassXiMinus    = 0.;
        Double_t lInvMassXiPlus     = 0.;
        Double_t lInvMassOmegaMinus = 0.;
        Double_t lInvMassOmegaPlus  = 0.;
        Bool_t lAssoXiMinus    = kFALSE;
        Bool_t lAssoXiPlus     = kFALSE;
        Bool_t lAssoOmegaMinus = kFALSE;
        Bool_t lAssoOmegaPlus  = kFALSE;
       
        Float_t  etaBach = 0.;
        Float_t  etaPos  = 0.;
        Float_t  etaNeg  = 0.;

        if ( fAnalysisType == "ESD" ) {		

             // - Load the cascade
             AliESDcascade *xiESD = lESDevent->GetCascade(iXi);
	     if (!xiESD) continue;
	
	     // - Connection to daughter tracks of the current cascade		
             UInt_t lIdxPosXi = (UInt_t) TMath::Abs( xiESD->GetPindex() );
             UInt_t lIdxNegXi = (UInt_t) TMath::Abs( xiESD->GetNindex() );
             UInt_t lBachIdx  = (UInt_t) TMath::Abs( xiESD->GetBindex() );
                   
             // - Rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
             if(lBachIdx == lIdxNegXi) {
                 AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
             }
             if(lBachIdx == lIdxPosXi) {
                 AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
             }
      
             // - Get the daughter tracks
	     AliESDtrack *pTrackXi = lESDevent->GetTrack( lIdxPosXi );
	     AliESDtrack *nTrackXi = lESDevent->GetTrack( lIdxNegXi );
	     AliESDtrack *bachTrackXi = lESDevent->GetTrack( lBachIdx  );
	     if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
	         Printf("ERROR: Could not retrieve one of the 3 daughter tracks of the cascade ...");
		 continue;
	     }
	
             // Get the number of TPC clusters
             lPosTPCClusters   = pTrackXi->GetTPCNcls();
             lNegTPCClusters   = nTrackXi->GetTPCNcls();
             lBachTPCClusters  = bachTrackXi->GetTPCNcls(); 
             // - Rejection of a poor quality tracks
             if(fkQualityCutTPCrefit){
                 // - Poor quality related to TPCrefit
                 ULong_t pStatus    = pTrackXi->GetStatus();
                 ULong_t nStatus    = nTrackXi->GetStatus();
                 ULong_t bachStatus = bachTrackXi->GetStatus();
                 if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
                 if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
                 if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
             } 
             if(fkQualityCutnTPCcls){
                // - Poor quality related to TPC clusters
                if(lPosTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Pos. track has less than 80 TPC clusters ... continue!"); continue; }
                if(lNegTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Neg. track has less than 80 TPC clusters ... continue!"); continue; }
                if(lBachTPCClusters < fMinnTPCcls) { AliWarning("Pb / Bach.   track has less than 80 TPC clusters ... continue!"); continue; }
             }

             etaPos  = pTrackXi->Eta();
             etaNeg  = nTrackXi->Eta();
             etaBach = bachTrackXi->Eta();
	
	     // - Info over reconstructed cascades
	     Double_t lV0quality = 0.;
	     if( bachTrackXi->Charge() < 0 ) {
                  //Calculate the effective mass of the Xi- candidate: Xi- hyp. (pdg code 3312
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
             lDcaXiDaughters          = xiESD->GetDcaXiDaughters();
             lDcaBachToPrimVertexXi   = TMath::Abs( bachTrackXi->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lMagneticField) );
             lXiCosineOfPointingAngle = xiESD->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
             xiESD->GetXYZcascade( lPosXi[0], lPosXi[1], lPosXi[2] );  
             lInvMassLambdaAsCascDghter = xiESD->GetEffMass();
             lDcaV0DaughtersXi          = xiESD->GetDcaV0Daughters();
             lV0CosineOfPointingAngleXi = xiESD->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
             lV0CosineOfPointingAngle   = xiESD->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2]); 
             xiESD->GetXYZ( lPosV0Xi[0], lPosV0Xi[1], lPosV0Xi[2] );
             lDcaV0ToPrimVertexXi       = xiESD->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
             lDcaPosToPrimVertexXi      = TMath::Abs( pTrackXi->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lMagneticField) );
             lDcaNegToPrimVertexXi      = TMath::Abs( nTrackXi->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lMagneticField) );
	     lChargeXi = xiESD->Charge();
	
             //------------------
	     // - PID Information

	     // - Combined VO-positive-daughter PID
	     AliPID pPidXi;    pPidXi.SetPriors( lPriorsGuessXi );
	     AliPID pPidOmega; pPidOmega.SetPriors( lPriorsGuessOmega );		
	     if( pTrackXi->IsOn(AliESDtrack::kESDpid) ){  
	          Double_t r[10] = {0.}; pTrackXi->GetESDpid(r);
		  pPidXi.SetProbabilities(r);
		  pPidOmega.SetProbabilities(r);		
		  // Check if the V0 positive track is a proton (case for Xi-)
		  Double_t pproton = pPidXi.GetProbability(AliPID::kProton);
		  if (pproton > pPidXi.GetProbability(AliPID::kElectron) &&
		      pproton > pPidXi.GetProbability(AliPID::kMuon)     &&
		      pproton > pPidXi.GetProbability(AliPID::kPion)     &&
		      pproton > pPidXi.GetProbability(AliPID::kKaon)        ) lIsPosInXiProton = kTRUE;
		  // Check if the V0 positive track is a pi+ (case for Xi+)
		  Double_t ppion = pPidXi.GetProbability(AliPID::kPion);
		  if (ppion > pPidXi.GetProbability(AliPID::kElectron) &&
		      ppion > pPidXi.GetProbability(AliPID::kMuon)     &&
		      ppion > pPidXi.GetProbability(AliPID::kKaon)     &&
		      ppion > pPidXi.GetProbability(AliPID::kProton)      ) lIsPosInXiPion = kTRUE;
		  // Check if the V0 positive track is a proton (case for Omega-)
		  pproton = 0.;
		  pproton = pPidOmega.GetProbability(AliPID::kProton);
		  if (pproton > pPidOmega.GetProbability(AliPID::kElectron) &&
		      pproton > pPidOmega.GetProbability(AliPID::kMuon)     &&
		      pproton > pPidOmega.GetProbability(AliPID::kPion)     &&
		      pproton > pPidOmega.GetProbability(AliPID::kKaon)       ) lIsPosInOmegaProton = kTRUE;
	 	  // Check if the V0 positive track is a pi+ (case for Omega+)
		  ppion = 0.;
		  ppion = pPidOmega.GetProbability(AliPID::kPion);
		  if (ppion > pPidOmega.GetProbability(AliPID::kElectron) &&
		      ppion > pPidOmega.GetProbability(AliPID::kMuon)     &&
		      ppion > pPidOmega.GetProbability(AliPID::kKaon)     &&
		      ppion > pPidOmega.GetProbability(AliPID::kProton)      ) lIsPosInOmegaPion = kTRUE;
	     }		
	     // - Combined VO-negative-daughter PID
	     AliPID nPidXi;    nPidXi.SetPriors( lPriorsGuessXi );
	     AliPID nPidOmega; nPidOmega.SetPriors( lPriorsGuessOmega );		
	     if( nTrackXi->IsOn(AliESDtrack::kESDpid) ) {  
     	          Double_t r[10] = {0.}; nTrackXi->GetESDpid(r);
	          nPidXi.SetProbabilities(r);
		  nPidOmega.SetProbabilities(r);
		  // Check if the V0 negative track is a pi- (case for Xi-)
		  Double_t ppion = nPidXi.GetProbability(AliPID::kPion);
		  if (ppion > nPidXi.GetProbability(AliPID::kElectron) &&
		      ppion > nPidXi.GetProbability(AliPID::kMuon)     &&
		      ppion > nPidXi.GetProbability(AliPID::kKaon)     &&
		      ppion > nPidXi.GetProbability(AliPID::kProton)      ) lIsNegInXiPion = kTRUE;
		  // Check if the V0 negative track is an anti-proton (case for Xi+)
		  Double_t pproton = nPidXi.GetProbability(AliPID::kProton);
		  if (pproton > nPidXi.GetProbability(AliPID::kElectron) &&
		      pproton > nPidXi.GetProbability(AliPID::kMuon)     &&
		      pproton > nPidXi.GetProbability(AliPID::kPion)     &&
		      pproton > nPidXi.GetProbability(AliPID::kKaon)       ) lIsNegInXiProton = kTRUE;
		  // Check if the V0 negative track is a pi- (case for Omega-)
		  ppion = 0.;
		  ppion = nPidOmega.GetProbability(AliPID::kPion);
		  if (ppion > nPidOmega.GetProbability(AliPID::kElectron) &&
		      ppion > nPidOmega.GetProbability(AliPID::kMuon)     &&
		      ppion > nPidOmega.GetProbability(AliPID::kKaon)     &&
		      ppion > nPidOmega.GetProbability(AliPID::kProton)     ) lIsNegInOmegaPion = kTRUE;
		  // Check if the V0 negative track is an anti-proton (case for Omega+)
		  pproton = 0.;
		  pproton = nPidOmega.GetProbability(AliPID::kProton);
		  if (pproton > nPidOmega.GetProbability(AliPID::kElectron) &&
		      pproton > nPidOmega.GetProbability(AliPID::kMuon)     &&
		      pproton > nPidOmega.GetProbability(AliPID::kPion)     &&
		      pproton > nPidOmega.GetProbability(AliPID::kKaon)       ) lIsNegInOmegaProton = kTRUE;
	     }
	     // - Combined bachelor PID
	     AliPID bachPidXi;	  bachPidXi.SetPriors( lPriorsGuessXi );
	     AliPID bachPidOmega; bachPidOmega.SetPriors( lPriorsGuessOmega );
	     if ( bachTrackXi->IsOn(AliESDtrack::kESDpid) ) {  
		   Double_t r[10] = {0.}; bachTrackXi->GetESDpid(r);
		   bachPidXi.SetProbabilities(r);
		   bachPidOmega.SetProbabilities(r);
		   // Check if the bachelor track is a pion
		   ppionBach = bachPidXi.GetProbability(AliPID::kPion);
		   if (ppionBach > bachPidXi.GetProbability(AliPID::kElectron) &&
		       ppionBach > bachPidXi.GetProbability(AliPID::kMuon)     &&
		       ppionBach > bachPidXi.GetProbability(AliPID::kKaon)     &&
		       ppionBach > bachPidXi.GetProbability(AliPID::kProton)     ) lIsBachelorPion = kTRUE;
		   // Check if the bachelor track is a kaon
		   pkaonBach = bachPidOmega.GetProbability(AliPID::kKaon);
		   if (pkaonBach > bachPidOmega.GetProbability(AliPID::kElectron) &&
		       pkaonBach > bachPidOmega.GetProbability(AliPID::kMuon)     &&
		       pkaonBach > bachPidOmega.GetProbability(AliPID::kPion)     &&
		       pkaonBach > bachPidOmega.GetProbability(AliPID::kProton)     ) lIsBachelorKaon = kTRUE;	
	     }
	     // - 4-sigma bands on Bethe-Bloch curve
             // Bachelor
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 4) lIsBachelorKaonForTPC = kTRUE;
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 4) lIsBachelorPionForTPC = kTRUE;
             // Negative V0 daughter
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 4) lIsNegPionForTPC   = kTRUE;
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 4) lIsNegProtonForTPC = kTRUE;
             // Positive V0 daughter
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 4) lIsPosPionForTPC   = kTRUE;
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 4) lIsPosProtonForTPC = kTRUE;
             /*        
             const AliExternalTrackParam *pInnerWallTrackXi    = pTrackXi    ->GetInnerParam(); // Do not use GetTPCInnerWall
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
             }*/
	     // - PID proba Vs Pt(Bach)
	     Int_t      lblBachForPID = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
	     TParticle* mcBachForPID  = lMCstack->Particle( lblBachForPID );
	     lmcPtBach = mcBachForPID->Pt();
	     // - MC perfect PID
	     if( mcBachForPID->GetPdgCode() == -211) lIsBachelorMCPiMinus = kTRUE;
	     if( mcBachForPID->GetPdgCode() ==  211) lIsBachelorMCPiPlus  = kTRUE;
	     if( mcBachForPID->GetPdgCode() == -321) lIsBachelorMCKMinus  = kTRUE;
	     if( mcBachForPID->GetPdgCode() ==  321) lIsBachelorMCKPlus   = kTRUE;
	
	
	     //---------------------------------------------------------------
	     // -  MC association (care : lots of "continue;" below this line)
	     if(fDebug > 5) cout<< "MC EventNumber: "<<lMCevent->Header()->GetEvent()<<" / MC event Number in Run : "<<lMCevent->Header()->GetEventNrInRun()<<endl;
	     // - Level of the V0 daughters
	     Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );  
	     Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );  		
	     TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
	     TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
	     // - Level of the Xi daughters	
	     Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ; 
	     Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();
             if( lblMotherPosV0Dghter != lblMotherNegV0Dghter) continue; // same mother
             if( lblMotherPosV0Dghter < 0 ) continue;                    // this particle is primary, no mother   
             if( lblMotherNegV0Dghter < 0 ) continue;                    // this particle is primary, no mother
	                                 	                                 // mothers = Lambda candidate ... a priori
	     TParticle* mcMotherPosV0Dghter = lMCstack->Particle( lblMotherPosV0Dghter );
	     TParticle* mcMotherNegV0Dghter = lMCstack->Particle( lblMotherNegV0Dghter );  // MN: redundant?? already checked that labels are the same...-->same part from stack
             Int_t lblBach = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
	     TParticle* mcBach = lMCstack->Particle( lblBach );	
	     // - Level of Xi candidate
	     Int_t lblGdMotherPosV0Dghter = mcMotherPosV0Dghter->GetFirstMother() ;
	     Int_t lblGdMotherNegV0Dghter = mcMotherNegV0Dghter->GetFirstMother() ;
             if( lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter ) continue;
             if( lblGdMotherPosV0Dghter < 0 ) continue;                  // primary lambda ...   
             if( lblGdMotherNegV0Dghter < 0 ) continue;                  // primary lambda ...   				
		                                                              // Gd mothers = Xi candidate ... a priori
             TParticle* mcGdMotherPosV0Dghter = lMCstack->Particle( lblGdMotherPosV0Dghter );
             TParticle* mcGdMotherNegV0Dghter = lMCstack->Particle( lblGdMotherNegV0Dghter );				
	     Int_t lblMotherBach = (Int_t) TMath::Abs( mcBach->GetFirstMother()  );  
             if( lblMotherBach != lblGdMotherPosV0Dghter ) continue; //same mother for bach and V0 daughters
	     TParticle* mcMotherBach = lMCstack->Particle( lblMotherBach );
            
             // - Check if cascade is primary
             if (!(lMCstack->IsPhysicalPrimary(lblMotherBach))) continue;  

	     // - Manage boolean for association
	     if     ( mcMotherBach          ->GetPdgCode() == 3312 &&
	              mcGdMotherPosV0Dghter ->GetPdgCode() == 3312 &&
	              mcGdMotherNegV0Dghter ->GetPdgCode() == 3312   ) {lAssoXiMinus = kTRUE;
                                                                        cascadeMass = 1.321;
                                                                         nAssoXiMinus++; }
	     else if( mcMotherBach           ->GetPdgCode() == -3312 &&
	              mcGdMotherPosV0Dghter  ->GetPdgCode() == -3312 &&
	              mcGdMotherNegV0Dghter  ->GetPdgCode() == -3312   ) {lAssoXiPlus = kTRUE;
                                                                          cascadeMass = 1.321;
                                                                          nAssoXiPlus++; }
	     else if( mcMotherBach           ->GetPdgCode() == 3334 &&
	              mcGdMotherPosV0Dghter  ->GetPdgCode() == 3334 &&
	              mcGdMotherNegV0Dghter  ->GetPdgCode() == 3334    ) {lAssoOmegaMinus = kTRUE;
                                                                          cascadeMass = 1.672;
                                                                          nAssoOmegaMinus++; }
	     else if( mcMotherBach           ->GetPdgCode() == -3334 &&
	              mcGdMotherPosV0Dghter  ->GetPdgCode() == -3334 &&
	              mcGdMotherNegV0Dghter  ->GetPdgCode() == -3334   ) {lAssoOmegaPlus = kTRUE;
                                                                             cascadeMass = 1.672;
                                                                             nAssoOmegaPlus++; }
	     // If a proper association  exists ...	
	     if(fDebug > 4){
	              cout<<"XiMinus    = "<<lAssoXiMinus   <<endl;
		      cout<<"XiPlus     = "<<lAssoXiPlus    <<endl;
		      cout<<"OmegaMinus = "<<lAssoOmegaMinus<<endl;
		      cout<<"OmegaPlus  = "<<lAssoOmegaPlus <<endl 
		          <<"----" 	      <<endl;	
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
             const AliAODcascade *xiAOD = lAODevent->GetCascade(iXi);
             if (!xiAOD) continue;

             // - Connection to daughter tracks of the current cascade
             AliAODTrack *pTrackXi    = dynamic_cast<AliAODTrack*>( xiAOD->GetDaughter(0) );
             AliAODTrack *nTrackXi    = dynamic_cast<AliAODTrack*>( xiAOD->GetDaughter(1) );
             AliAODTrack *bachTrackXi = dynamic_cast<AliAODTrack*>( xiAOD->GetDecayVertexXi()->GetDaughter(0) );
             if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
                  AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
                  continue;
             }
             UInt_t lIdxPosXi  = (UInt_t) TMath::Abs( pTrackXi->GetID() );
             UInt_t lIdxNegXi  = (UInt_t) TMath::Abs( nTrackXi->GetID() );
             UInt_t lBachIdx   = (UInt_t) TMath::Abs( bachTrackXi->GetID() );

             // - Rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
             if(lBachIdx == lIdxNegXi) {
                AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
             }
             if(lBachIdx == lIdxPosXi) {
                AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
             }
             lPosTPCClusters   = pTrackXi->GetTPCNcls();
             lNegTPCClusters   = nTrackXi->GetTPCNcls();
             lBachTPCClusters  = bachTrackXi->GetTPCNcls();

             // - Rejection of a poor quality tracks
             if (fkQualityCutTPCrefit) {
                 // - Poor quality related to TPCrefit
                 if (!(pTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
                 if (!(nTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
                 if (!(bachTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
             }
             if (fkQualityCutnTPCcls) {
                 // - Poor quality related to TPC clusters
                 if(lPosTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Pos. track has less than 80 TPC clusters ... continue!"); continue; }
                 if(lNegTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Neg. track has less than 80 TPC clusters ... continue!"); continue; }
                 if(lBachTPCClusters < fMinnTPCcls) { AliWarning("Pb / Bach.   track has less than 80 TPC clusters ... continue!"); continue; }
             }

             etaPos  = pTrackXi->Eta();
             etaNeg  = nTrackXi->Eta();
             etaBach = bachTrackXi->Eta();

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

             //------------------
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
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 4) lIsBachelorKaonForTPC = kTRUE;
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 4) lIsBachelorPionForTPC = kTRUE;
             // Negative V0 daughter
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 4) lIsNegPionForTPC   = kTRUE;
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 4) lIsNegProtonForTPC = kTRUE;
             // Positive V0 daughter
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 4) lIsPosPionForTPC   = kTRUE;
             if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 4) lIsPosProtonForTPC = kTRUE;
             /*
             const AliExternalTrackParam *pInnerWallTrackXi    = pTrackXi    ->GetInnerParam(); // Do not use GetTPCInnerWall
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
             }*/

             // - PID proba Vs Pt(Bach)
             Int_t lblBachForPID = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
             AliAODMCParticle* mcBachForPID   = (AliAODMCParticle*) arrayMC->At( lblBachForPID );
             lmcPtBach = mcBachForPID->Pt();

             // - MC perfect PID
             if( mcBachForPID->PdgCode() == -211) lIsBachelorMCPiMinus = kTRUE;
             if( mcBachForPID->PdgCode() ==  211) lIsBachelorMCPiPlus  = kTRUE;
             if( mcBachForPID->PdgCode() == -321) lIsBachelorMCKMinus  = kTRUE;
             if( mcBachForPID->PdgCode() ==  321) lIsBachelorMCKPlus   = kTRUE;

             //--------------------------------------------------------------
             // - MC association (care : lots of "continue;" below this line)
             if(fDebug > 5) cout<<"MC EventNumber : "<<lMCevent->Header()->GetEvent()<<" / MC event Number in Run : "<<lMCevent->Header()->GetEventNrInRun()<<endl;
             // - Level of the V0 daughters
             Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );
             Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );
             AliAODMCParticle* mcPosV0Dghter = (AliAODMCParticle*) arrayMC->At( lblPosV0Dghter );
             AliAODMCParticle* mcNegV0Dghter = (AliAODMCParticle*) arrayMC->At( lblNegV0Dghter );
             // - Level of the Xi daughters
             Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetMother();   
             Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetMother();
             if( lblMotherPosV0Dghter != lblMotherNegV0Dghter) continue; // same mother
             if( lblMotherPosV0Dghter < 0 ) continue;                    // this particle is primary, no mother
             if( lblMotherNegV0Dghter < 0 ) continue;                    // this particle is primary, no mother
                                                                         // mothers = Lambda candidate ... a priori
             AliAODMCParticle* mcMotherPosV0Dghter = (AliAODMCParticle*) arrayMC->At( lblMotherPosV0Dghter );
             AliAODMCParticle* mcMotherNegV0Dghter = (AliAODMCParticle*) arrayMC->At( lblMotherNegV0Dghter );  
             Int_t      lblBach  = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
             AliAODMCParticle* mcBach   = (AliAODMCParticle*) arrayMC->At( lblBach );
             // - Level of Xi candidate
             Int_t lblGdMotherPosV0Dghter =   mcMotherPosV0Dghter->GetMother() ;
             Int_t lblGdMotherNegV0Dghter =   mcMotherNegV0Dghter->GetMother() ;
             if( lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter ) continue;
             if( lblGdMotherPosV0Dghter < 0 ) continue;                    // primary lambda ...
             if( lblGdMotherNegV0Dghter < 0 ) continue;                    // primary lambda ...
                                                                           // Gd mothers = Xi candidate ... a priori
             AliAODMCParticle* mcGdMotherPosV0Dghter = (AliAODMCParticle*) arrayMC->At( lblGdMotherPosV0Dghter );
             AliAODMCParticle* mcGdMotherNegV0Dghter = (AliAODMCParticle*) arrayMC->At( lblGdMotherNegV0Dghter );
             Int_t lblMotherBach = (Int_t) TMath::Abs( mcBach->GetMother() );
             if( lblMotherBach != lblGdMotherPosV0Dghter ) continue; //same mother for bach and V0 daughters
             AliAODMCParticle* mcMotherBach = (AliAODMCParticle*) arrayMC->At( lblMotherBach );

             // - Check if cascade is primary
             if (!(mcMotherBach->IsPhysicalPrimary())) continue;

             // - Manage boolean for association
             if     ( mcMotherBach          ->GetPdgCode() == 3312 &&
                      mcGdMotherPosV0Dghter ->GetPdgCode() == 3312 &&
                      mcGdMotherNegV0Dghter ->GetPdgCode() == 3312   ) {lAssoXiMinus = kTRUE;
                                                                        cascadeMass = 1.321;
                                                                        nAssoXiMinus++; }
             else if( mcMotherBach           ->GetPdgCode() == -3312 &&
                      mcGdMotherPosV0Dghter  ->GetPdgCode() == -3312 &&
                      mcGdMotherNegV0Dghter  ->GetPdgCode() == -3312   ) {lAssoXiPlus = kTRUE;
                                                                          cascadeMass = 1.321;
                                                                          nAssoXiPlus++; }
             else if( mcMotherBach           ->GetPdgCode() == 3334 &&
                      mcGdMotherPosV0Dghter  ->GetPdgCode() == 3334 &&
                      mcGdMotherNegV0Dghter  ->GetPdgCode() == 3334   ) {lAssoOmegaMinus = kTRUE;
                                                                         cascadeMass = 1.672;
                                                                         nAssoOmegaMinus++; }
             else if( mcMotherBach           ->GetPdgCode() == -3334 &&
                      mcGdMotherPosV0Dghter  ->GetPdgCode() == -3334 &&
                      mcGdMotherNegV0Dghter  ->GetPdgCode() == -3334   ) {lAssoOmegaPlus = kTRUE;
                                                                          cascadeMass = 1.672;
                                                                          nAssoOmegaPlus++; }

             lmcPt              = mcMotherBach->Pt();
             lmcRapCasc         = 0.5*TMath::Log( (mcMotherBach->E() + mcMotherBach->Pz()) / (mcMotherBach->E() - mcMotherBach->Pz() +1.e-13) );
             lmcEta             = mcMotherBach->Eta();
             Float_t decayCascX = mcBach->Xv();
             Float_t decayCascY = mcBach->Yv();
             lmcTransvRadius    = TMath::Sqrt(decayCascX*decayCascX+decayCascY*decayCascY); // decay point of Xi, = the production vertex of Bachelor ...
             TVector3 lmcTVect3Mom( mcMotherBach->Px(), mcMotherBach->Py(), mcMotherBach->Pz() );
             Double_t xiMomX    = xiAOD->MomXiX();
             Double_t xiMomY    = xiAOD->MomXiY();
             Double_t xiMomZ    = xiAOD->MomXiZ();
             lrecoPt            = TMath::Sqrt( xiMomX*xiMomX   + xiMomY*xiMomY ); 
             lrecoTransvRadius  = TMath::Sqrt( xiAOD->DecayVertexXiX() * xiAOD->DecayVertexXiX() + xiAOD->DecayVertexXiY() * xiAOD->DecayVertexXiY() );
             TVector3 lrecoTVect3Mom( xiMomX, xiMomY, xiMomZ );
             lDeltaPhiMcReco    = lmcTVect3Mom.DeltaPhi( lrecoTVect3Mom ) * 180.0/TMath::Pi();
             lmcPtPosV0Dghter   = mcPosV0Dghter->Pt() ;
             lmcPtNegV0Dghter   = mcNegV0Dghter->Pt();
             lrecoP             = TMath::Sqrt( xiMomX*xiMomX   + xiMomY*xiMomY   + xiMomZ*xiMomZ );;
             Double_t lV0momX   = xiAOD->MomV0X();
             Double_t lV0momY   = xiAOD->MomV0Y();
             Double_t lV0momZ   = xiAOD->MomV0Z();
             lV0mom             = TMath::Sqrt(TMath::Power(lV0momX,2)+TMath::Power(lV0momY,2)+TMath::Power(lV0momZ,2));
             Double_t lBachMomX = xiAOD->MomBachX();
             Double_t lBachMomY = xiAOD->MomBachY();
             lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
             Double_t lV0NMomX = xiAOD->MomNegX();
             Double_t lV0NMomY = xiAOD->MomNegY();
             Double_t lV0PMomX = xiAOD->MomPosX();
             Double_t lV0PMomY = xiAOD->MomPosY();
             lnTrackTransvMom = TMath::Sqrt( lV0NMomX*lV0NMomX   + lV0NMomY*lV0NMomY );
             lpTrackTransvMom = TMath::Sqrt( lV0PMomX*lV0PMomX   + lV0PMomY*lV0PMomY );
            
        }

        lXiRadius   = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
        lV0RadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] ); 

        // - Cut on pt of the three daughter tracks
        if (lBachTransvMom<fMinPtCutOnDaughterTracks) continue;
        if (lpTrackTransvMom<fMinPtCutOnDaughterTracks) continue;
        if (lnTrackTransvMom<fMinPtCutOnDaughterTracks) continue;
       
        // - Cut on pseudorapidity of the three daughter tracks
        if (TMath::Abs(etaBach)>fEtaCutOnDaughterTracks) continue;
        if (TMath::Abs(etaPos)>fEtaCutOnDaughterTracks) continue;
        if (TMath::Abs(etaNeg)>fEtaCutOnDaughterTracks) continue;
       
        // - Extra-selection for cascade candidates
        if (fkExtraSelections) {
                if (lDcaXiDaughters > 0.3) continue;              // in AliCascadeVertexer
                if (lXiCosineOfPointingAngle < 0.999 ) continue;  // in AliCascadeVertexer
                if (lDcaV0ToPrimVertexXi < 0.05) continue;        // in AliCascadeVertexer
                if (lDcaBachToPrimVertexXi < 0.03) continue;      // in AliCascadeVertexer
                if (lDcaV0DaughtersXi > 1.) continue;             // in AliV0vertexer
                if (lV0CosineOfPointingAngleXi < 0.998) continue; // in AliV0vertexer
                if (lDcaPosToPrimVertexXi < 0.1) continue;        // in AliV0vertexer
                if (lDcaNegToPrimVertexXi < 0.1) continue;        // in AliV0vertexer
                if(lXiRadius < .9) continue;                      // in AliCascadeVertexer
                if(lV0RadiusXi < 0.9) continue;                   // in AliV0vertexer
        }

        //-------------------------
        // - Fill combined PID TH1s
        if( lChargeXi < 0 && lIsBachelorPion )    fHistMassWithCombPIDXiMinus    ->Fill( lInvMassXiMinus    );
        if( lChargeXi > 0 && lIsBachelorPion )    fHistMassWithCombPIDXiPlus     ->Fill( lInvMassXiPlus     );
        if( lChargeXi < 0 && lIsBachelorKaon )    fHistMassWithCombPIDOmegaMinus ->Fill( lInvMassOmegaMinus );
        if( lChargeXi > 0 && lIsBachelorKaon )    fHistMassWithCombPIDOmegaPlus  ->Fill( lInvMassOmegaPlus  );
        if( lChargeXi < 0 )   fHistMassXiMinus    ->Fill( lInvMassXiMinus );
        if( lChargeXi > 0 )   fHistMassXiPlus     ->Fill( lInvMassXiPlus );
        if( lChargeXi < 0 )   fHistMassOmegaMinus ->Fill( lInvMassOmegaMinus );
        if( lChargeXi > 0 )   fHistMassOmegaPlus  ->Fill( lInvMassOmegaPlus );
        if(lIsBachelorPion)   f2dHistPIDprobaPionVsMCPtBach->Fill( lmcPtBach, ppionBach );
        if(lIsBachelorKaon)   f2dHistPIDprobaKaonVsMCPtBach->Fill( lmcPtBach, pkaonBach );
        if( lChargeXi < 0 && lIsBachelorMCPiMinus )    fHistMassWithMcPIDXiMinus     ->Fill( lInvMassXiMinus );
        if( lChargeXi > 0 && lIsBachelorMCPiPlus  )    fHistMassWithMcPIDXiPlus      ->Fill( lInvMassXiPlus );
        if( lChargeXi < 0 && lIsBachelorMCKMinus  )    fHistMassWithMcPIDOmegaMinus  ->Fill( lInvMassOmegaMinus );
        if( lChargeXi > 0 && lIsBachelorMCKPlus   )    fHistMassWithMcPIDOmegaPlus   ->Fill( lInvMassOmegaPlus );


        // - No association, skip the rest of the code
        if(!lAssoXiMinus && !lAssoXiPlus && !lAssoOmegaMinus && !lAssoOmegaPlus) continue; 

        //--------------
        // - Proper time         
        // For cascade (reconstructed)   
        Double_t lctau = TMath::Sqrt(TMath::Power((lPosXi[0]-lBestPrimaryVtxPos[0]),2)+TMath::Power((lPosXi[1]-lBestPrimaryVtxPos[1]),2)+TMath::Power((lPosXi[2]-lBestPrimaryVtxPos[2]),2));
        if (lrecoP!=0) lctau = lctau*cascadeMass/lrecoP;   
        else lctau = -1.;
        // For Lambda (reconstructed)
        Float_t lambdaMass = 1.115683; // PDG mass 
        Float_t distV0Xi = TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2)+TMath::Power((lPosV0Xi[2]-lPosXi[2]),2)); 
        Float_t lctauV0 = -1.;
        if (lV0mom!=0) lctauV0 = distV0Xi*lambdaMass/lV0mom; 
        // Distance
        Float_t distTV0Xi = TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2));
       
        //------------------------------------------------------------
        // - Fill histos for the cascade candidates associated with MC
        if( lChargeXi < 0 && lAssoXiMinus){	
                fHistAsMCMassXiMinus	      ->Fill( lInvMassXiMinus  );
		if(lIsBachelorPion)	f2dHistAsMCandCombPIDGenPtVsGenYXiMinus->Fill( lmcPt, lmcRapCasc );
		f2dHistAsMCGenPtVsGenYXiMinus ->Fill( lmcPt, lmcRapCasc);
		fHistAsMCGenEtaXiMinus        ->Fill( lmcEta           );
		f2dHistAsMCResPtXiMinus       ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResRXiMinus        ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
                f2dHistAsMCResPhiXiMinus      ->Fill( lmcPt, lDeltaPhiMcReco );
                f2dHistAsMCptProtonMCptXiMinus->Fill(lmcPt,lmcPtPosV0Dghter);
                fHistV0CosineOfPointingAnglevsPtXi->Fill(lmcPt,lV0CosineOfPointingAngle);
        }	
        else if( lChargeXi > 0 && lAssoXiPlus){	
		fHistAsMCMassXiPlus	      ->Fill( lInvMassXiPlus   );
		if(lIsBachelorPion)	f2dHistAsMCandCombPIDGenPtVsGenYXiPlus->Fill( lmcPt, lmcRapCasc );
		f2dHistAsMCGenPtVsGenYXiPlus  ->Fill( lmcPt, lmcRapCasc);
		fHistAsMCGenEtaXiPlus         ->Fill( lmcEta           );
		f2dHistAsMCResPtXiPlus        ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResRXiPlus         ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
                f2dHistAsMCResPhiXiPlus       ->Fill( lmcPt, lDeltaPhiMcReco );
                f2dHistAsMCptAntiprotonMCptXiPlus->Fill(lmcPt,lmcPtNegV0Dghter);
                fHistV0CosineOfPointingAnglevsPtXi->Fill(lmcPt,lV0CosineOfPointingAngle);
        }
        else if( lChargeXi < 0 && lAssoOmegaMinus){	
		fHistAsMCMassOmegaMinus          ->Fill( lInvMassOmegaMinus );
		if(lIsBachelorKaon)	f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus->Fill( lmcPt, lmcRapCasc );
		f2dHistAsMCGenPtVsGenYOmegaMinus ->Fill( lmcPt, lmcRapCasc  );
		fHistAsMCGenEtaOmegaMinus        ->Fill( lmcEta             );
		f2dHistAsMCResPtOmegaMinus       ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResROmegaMinus        ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
                f2dHistAsMCResPhiOmegaMinus      ->Fill( lmcPt, lDeltaPhiMcReco );
                f2dHistAsMCptProtonMCptOmegaMinus->Fill(lmcPt,lmcPtPosV0Dghter);
                fHistV0CosineOfPointingAnglevsPtOmega->Fill(lmcPt,lV0CosineOfPointingAngle);
        }	
        else if( lChargeXi > 0 && lAssoOmegaPlus){	
		fHistAsMCMassOmegaPlus           ->Fill( lInvMassOmegaPlus );
		if(lIsBachelorKaon)	f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus->Fill( lmcPt, lmcRapCasc );
		f2dHistAsMCGenPtVsGenYOmegaPlus  ->Fill( lmcPt, lmcRapCasc   );
		fHistAsMCGenEtaOmegaPlus         ->Fill( lmcEta            );
		f2dHistAsMCResPtOmegaPlus        ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResROmegaPlus         ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
                f2dHistAsMCResPhiOmegaPlus       ->Fill( lmcPt, lDeltaPhiMcReco );
                f2dHistAsMCptAntiprotonMCptOmegaPlus->Fill(lmcPt,lmcPtNegV0Dghter);
                fHistV0CosineOfPointingAnglevsPtOmega->Fill(lmcPt,lV0CosineOfPointingAngle);
        }
        fHistV0toXiCosineOfPointingAngle->Fill(lV0CosineOfPointingAngleXi);

        //------------------         
        // - Fill containers
       
        // - Filling the AliCFContainer (optimisation of topological selections + systematics)
        Double_t lContainerCutVars[19] = {0.0};
        lContainerCutVars[0]  = lDcaXiDaughters;
        lContainerCutVars[1]  = lDcaBachToPrimVertexXi;
        lContainerCutVars[2]  = lXiCosineOfPointingAngle;
        lContainerCutVars[3]  = lXiRadius;
        lContainerCutVars[4]  = lInvMassLambdaAsCascDghter;
        lContainerCutVars[5]  = lDcaV0DaughtersXi;
        lContainerCutVars[6]  = lV0CosineOfPointingAngleXi;
        lContainerCutVars[7]  = lV0RadiusXi;
        lContainerCutVars[8]  = lDcaV0ToPrimVertexXi;	
        lContainerCutVars[9]  = lDcaPosToPrimVertexXi;
        lContainerCutVars[10] = lDcaNegToPrimVertexXi;
        lContainerCutVars[13] = lmcPt;
        lContainerCutVars[16] = lctau;
        lContainerCutVars[17] = lctauV0;
        lContainerCutVars[18] = distTV0Xi;
        // All cases should be covered below
        if( lChargeXi < 0 && lAssoXiMinus    ) {
                     lContainerCutVars[11] = lInvMassXiMinus;
                     lContainerCutVars[12] = lInvMassOmegaMinus;//1.63;
                     lContainerCutVars[14] = lmcRapCasc;
                     lContainerCutVars[15] = -1.;
                     if ( lIsBachelorPionForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )    
                       fCFContAsCascadeCuts->Fill(lContainerCutVars,0); // for Xi-
        }
        if( lChargeXi > 0 && lAssoXiPlus    ) {
                     lContainerCutVars[11] = lInvMassXiPlus;
                     lContainerCutVars[12] = lInvMassOmegaPlus;//1.26;
                     lContainerCutVars[14] = lmcRapCasc;
                     lContainerCutVars[15] = -1.; 
                     if ( lIsBachelorPionForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )    
                       fCFContAsCascadeCuts->Fill(lContainerCutVars,1); // for Xi+
        }
        if( lChargeXi < 0 && lAssoOmegaMinus ) {
                     lContainerCutVars[11] = lInvMassXiMinus;//1.63;
                     lContainerCutVars[12] = lInvMassOmegaMinus;
                     lContainerCutVars[14] = -1.;
                     lContainerCutVars[15] = lmcRapCasc;
                     if ( lIsBachelorKaonForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )    
                       fCFContAsCascadeCuts->Fill(lContainerCutVars,2); // for Omega-
        }
        if( lChargeXi > 0 && lAssoOmegaPlus  ) {
                     lContainerCutVars[11] = lInvMassXiPlus;//1.26;
                     lContainerCutVars[12] = lInvMassOmegaPlus;
                     lContainerCutVars[14] = -1.;
                     lContainerCutVars[15] = lmcRapCasc;
                     if ( lIsBachelorKaonForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )    
                       fCFContAsCascadeCuts->Fill(lContainerCutVars,3); // for Omega+
        }
        
        // - Filling the AliCFContainers related to PID
        Double_t lContainerPIDVars[3] = {0.0};
	
        // Xi Minus		
        if( lChargeXi < 0 && lAssoXiMinus ) {
	        lContainerPIDVars[0] = lmcPt;
	        lContainerPIDVars[1] = lInvMassXiMinus;
	        lContainerPIDVars[2] = lmcRapCasc;
	        // No PID
	        fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 0); // No PID
	        // TPC PID
	        if( lIsBachelorPionForTPC ) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
	        if( lIsBachelorPionForTPC && lIsPosProtonForTPC ) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks	
	        if( lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	        // Combined PID
	        if( lIsBachelorPion ) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
	        if( lIsBachelorPion && lIsPosInXiProton ) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
	        if( lIsBachelorPion && lIsPosInXiProton && lIsNegInXiPion ) fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
        }	
        lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; 
	
        // Xi Plus		
        if( lChargeXi > 0 && lAssoXiPlus ) {
               lContainerPIDVars[0] = lmcPt;
               lContainerPIDVars[1] = lInvMassXiPlus;
     	       lContainerPIDVars[2] = lmcRapCasc;
	       // No PID
	       fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 0); // No PID
	       // TPC PID
	       if( lIsBachelorPionForTPC ) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
	       if( lIsBachelorPionForTPC && lIsNegProtonForTPC ) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
	       if( lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC ) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	       // Combined PID
	       if( lIsBachelorPion ) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
	       if( lIsBachelorPion && lIsNegInXiProton ) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
	       if( lIsBachelorPion && lIsNegInXiProton && lIsPosInXiPion ) fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
        }	
        lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; 
	
        // Omega Minus		
        if( lChargeXi < 0 && lAssoOmegaMinus ) {
	       lContainerPIDVars[0] = lmcPt;
	       lContainerPIDVars[1] = lInvMassOmegaMinus;
	       lContainerPIDVars[2] = lmcRapCasc;		
	       // No PID
	       fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 0); // No PID
	       // TPC PID
	       if( lIsBachelorKaonForTPC ) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
	       if( lIsBachelorKaonForTPC && lIsPosProtonForTPC ) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
	       if( lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	       // Combined PID
	       if( lIsBachelorKaon ) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
	       if( lIsBachelorKaon && lIsPosInOmegaProton ) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
	       if( lIsBachelorKaon && lIsPosInOmegaProton && lIsNegInOmegaPion ) fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
        }	
        lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; 
	
        // Omega Plus		
        if( lChargeXi > 0 && lAssoOmegaPlus) {
	       lContainerPIDVars[0] = lmcPt;
	       lContainerPIDVars[1] = lInvMassOmegaPlus;
	       lContainerPIDVars[2] = lmcRapCasc;		
	       // No PID
	       fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 0); // No PID
	       // TPC PID
	       if( lIsBachelorKaonForTPC ) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		     if( lIsBachelorKaonForTPC && lIsNegProtonForTPC ) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		     if( lIsBachelorKaonForTPC && lIsNegProtonForTPC && lIsPosPionForTPC ) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		     // Combined PID
		     if( lIsBachelorKaon ) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		     if( lIsBachelorKaon && lIsNegInOmegaProton ) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		     if( lIsBachelorKaon && lIsNegInOmegaProton && lIsPosInOmegaPion ) fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
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
void AliAnalysisTaskCheckPerformanceCascadepp276::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
	
 /* TList *cRetrievedList = 0x0;
  cRetrievedList = (TList*)GetOutputData(1);
  if(!cRetrievedList) {
	Printf("ERROR - AliAnalysisTaskCheckPerformanceCascadepp276 : ouput data container list not available\n");
	return;
  }	
	
  fHistTrackMultiplicityBeforeAnySel = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistTrackMultiplicityBeforeAnySel")  );
  if (!fHistTrackMultiplicityBeforeAnySel) {
    Printf("ERROR - AliAnalysisTaskCheckPerformanceCascadepp276 : fHistTrackMultiplicityBeforeAnySel not available");
    return;
  }
  
   
  TCanvas *canCheckPerformanceCascade = new TCanvas("AliAnalysisTaskCheckPerformanceCascadepp276","Multiplicity",10,10,510,510);
  canCheckPerformanceCascade->cd(1)->SetLogy();

  fHistTrackMultiplicityBeforeAnySel->SetMarkerStyle(22);
  fHistTrackMultiplicityBeforeAnySel->DrawCopy("E");
 */
}
