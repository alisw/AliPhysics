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
//                        - cut on TPC clusters as a parameter
//                        - cut on min pt of daughter tracks added (parameter+control histos)
//                        - cut on pseudorapidity for daughter tracks as a parameter (+control histos for Xi-)  
//               December2014
//                - de-commented the v0 and cascade revertexers
//                - added methods to define the set of cuts used in the V0 and cascade re-vertexers
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

#include "AliCentrality.h"

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
#include "AliAnalysisTaskCheckPerformanceCascadePbPb.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskCheckPerformanceCascadePbPb)



     //_____Dummy constructor________________________________________________________________
AliAnalysisTaskCheckPerformanceCascadePbPb::AliAnalysisTaskCheckPerformanceCascadePbPb() 
: AliAnalysisTaskSE(), // <- take care to AliAnalysisTask( empty )
  fAnalysisType("ESD"), fESDtrackCuts(0), /*fPaveTextBookKeeping(0),*/

    fPIDResponse                   (0),
    fkRerunV0CascVertexers         (0),
    fkQualityCutZprimVtxPos        (kTRUE),
    fkRejectEventPileUp            (kTRUE),
    fkQualityCutNoTPConlyPrimVtx   (kTRUE),
    fkQualityCutTPCrefit           (kTRUE),
    fkQualityCutnTPCcls            (kTRUE),
    fMinnTPCcls                    (0), 
    fkExtraSelections              (0),
    fCentrLowLim(0),    fCentrUpLim(0), fCentrEstimator(0), fkUseCleaning(0),
    fVtxRange                      (0),
    fApplyAccCut                   (0),
    fMinPtCutOnDaughterTracks      (0),
    fEtaCutOnDaughterTracks      (0),
    
	// - Cascade part initialisation
    fListHistCascade(0),

    // Events in centrality bins
    fHistEvtsInCentralityBinsvsNtracks(0),
    fHistBestVtxX(0),
    fHistBestVtxY(0),
    fHistBestVtxZ(0),
 
    // Cascade multiplicity histos

    fHistnXiPlusPerEvTot(0),
    fHistnXiMinusPerEvTot(0),
    fHistnOmegaPlusPerEvTot(0),
    fHistnOmegaMinusPerEvTot(0),

    fHistnXiPlusPerEv(0),  
    fHistnXiMinusPerEv(0),
    fHistnOmegaPlusPerEv(0),
    fHistnOmegaMinusPerEv(0),

    fHistnAssoXiMinus(0),
    fHistnAssoXiPlus(0),
    fHistnAssoOmegaMinus(0),
    fHistnAssoOmegaPlus(0),


    fHistMCTrackMultiplicity(0), 
       // - Resolution of the multiplicity estimator
    f2dHistRecoMultVsMCMult(0),

    fHistEtaGenProton(0),
    fHistEtaGenAntiProton(0),
      
   // Xi-
   fHistEtaGenCascXiMinus(0),
   f3dHistGenPtVsGenYvsCentXiMinusNat(0),
   f3dHistGenPtVsGenYvsNtracksXiMinusNat(0),
   f3dHistGenPtVsGenYvsCentXiMinusInj(0),
   f3dHistGenPtVsGenYvsNtracksXiMinusInj(0),
   f3dHistGenPtVsGenctauvsCentXiMinusNat(0),
   f3dHistGenPtVsGenctauvsCentXiMinusInj(0),
   
    fHistThetaGenCascXiMinusNat(0),
    fHistThetaGenCascXiMinusInj(0), 
    f2dHistGenPtVsGenYFdblXiMinus(0),
    
    fHistThetaLambdaXiMinus(0), 
    fHistThetaBachXiMinus(0),
    
    fHistThetaMesDghterXiMinus(0), 
    fHistThetaBarDghterXiMinus(0),
    
    fHistPtBachXiMinus(0),
    fHistPtMesDghterXiMinus(0),
    fHistPtBarDghterXiMinus(0),
   
    fHistPtRecBachXiMinus(0),
    fHistPtRecMesDghterXiMinus(0),
    fHistPtRecBarDghterXiMinus(0),

   // Xi+
   fHistEtaGenCascXiPlus(0),
   f3dHistGenPtVsGenYvsCentXiPlusNat(0),
   f3dHistGenPtVsGenYvsNtracksXiPlusNat(0),
   f3dHistGenPtVsGenYvsCentXiPlusInj(0),
   f3dHistGenPtVsGenYvsNtracksXiPlusInj(0),
   f3dHistGenPtVsGenctauvsCentXiPlusNat(0),
   f3dHistGenPtVsGenctauvsCentXiPlusInj(0),
   
    fHistThetaGenCascXiPlusNat(0), 
    fHistThetaGenCascXiPlusInj(0),
    f2dHistGenPtVsGenYFdblXiPlus(0),
    
    fHistThetaLambdaXiPlus(0), 
    fHistThetaBachXiPlus(0),
    
    fHistThetaMesDghterXiPlus(0), 
    fHistThetaBarDghterXiPlus(0),
    
    fHistPtBachXiPlus(0),
    fHistPtMesDghterXiPlus(0),
    fHistPtBarDghterXiPlus(0),
  
   // Omega-
   fHistEtaGenCascOmegaMinus(0),
   f3dHistGenPtVsGenYvsCentOmegaMinusNat(0),
   f3dHistGenPtVsGenYvsNtracksOmegaMinusNat(0),
   f3dHistGenPtVsGenYvsCentOmegaMinusInj(0),
   f3dHistGenPtVsGenYvsNtracksOmegaMinusInj(0),
   f3dHistGenPtVsGenctauvsCentOmegaMinusNat(0),
   f3dHistGenPtVsGenctauvsCentOmegaMinusInj(0),
   
    fHistThetaGenCascOmegaMinusNat(0),
    fHistThetaGenCascOmegaMinusInj(0),
    f2dHistGenPtVsGenYFdblOmegaMinus(0),
    
    fHistThetaLambdaOmegaMinus(0), 
    fHistThetaBachOmegaMinus(0),
    
    fHistThetaMesDghterOmegaMinus(0), 
    fHistThetaBarDghterOmegaMinus(0),
    
    fHistPtBachOmegaMinus(0),
    fHistPtMesDghterOmegaMinus(0),
    fHistPtBarDghterOmegaMinus(0),
   
   // Omega+
   fHistEtaGenCascOmegaPlus(0),
   f3dHistGenPtVsGenYvsCentOmegaPlusNat(0),
   f3dHistGenPtVsGenYvsNtracksOmegaPlusNat(0),
   f3dHistGenPtVsGenYvsCentOmegaPlusInj(0),
   f3dHistGenPtVsGenYvsNtracksOmegaPlusInj(0),
   f3dHistGenPtVsGenctauvsCentOmegaPlusNat(0),
   f3dHistGenPtVsGenctauvsCentOmegaPlusInj(0),
   
    fHistThetaGenCascOmegaPlusNat(0),
    fHistThetaGenCascOmegaPlusInj(0),
    f2dHistGenPtVsGenYFdblOmegaPlus(0),
    
    fHistThetaLambdaOmegaPlus(0), 
    fHistThetaBachOmegaPlus(0),
    
    fHistThetaMesDghterOmegaPlus(0), 
    fHistThetaBarDghterOmegaPlus(0),
    
    fHistPtBachOmegaPlus(0),
    fHistPtMesDghterOmegaPlus(0),
    fHistPtBarDghterOmegaPlus(0),

// Part 2 - Association to MC
	
    fHistMassXiMinus(0),
    fHistMassXiPlus(0),
    fHistMassOmegaMinus(0),
    fHistMassOmegaPlus(0),
    
	// - Effective mass histos with combined PID
    fHistMassWithCombPIDXiMinus(0), fHistMassWithCombPIDXiPlus(0),
    fHistMassWithCombPIDOmegaMinus(0), fHistMassWithCombPIDOmegaPlus(0),
    	
	// - PID Probability versus MC Pt(bachelor track)
    f2dHistPIDprobaKaonVsMCPtBach(0), f2dHistPIDprobaPionVsMCPtBach(0),
    
    	// - Effective mass histos with perfect MC PID on the bachelor
    fHistMassWithMcPIDXiMinus(0), fHistMassWithMcPIDXiPlus(0),
    fHistMassWithMcPIDOmegaMinus(0), fHistMassWithMcPIDOmegaPlus(0),
	
    
	// - Effective mass histos for the cascade candidates associated with MC
    fHistAsMCMassXiMinus(0),		
    fHistAsMCMassXiPlus(0),		
    fHistAsMCMassOmegaMinus(0),
    fHistAsMCMassOmegaPlus(0),
    
	// - Generated Pt Vs generated y, for the cascade candidates associated with MC + Info Comb. PID
    f2dHistAsMCandCombPIDGenPtVsGenYXiMinus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYXiPlus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus(0),
	   
    	// - Generated Pt Vs generated y, for the cascade candidates associated with MC
    f2dHistAsMCGenPtVsGenYXiMinus(0),
    f2dHistAsMCGenPtVsGenYXiPlus(0),
    f2dHistAsMCGenPtVsGenYOmegaMinus(0),
    f2dHistAsMCGenPtVsGenYOmegaPlus(0),
    
    	// - Generated Eta of the the cascade candidates associated with MC
    fHistAsMCGenEtaXiMinus(0),
    fHistAsMCGenEtaXiPlus(0),
    fHistAsMCGenEtaOmegaMinus(0),
    fHistAsMCGenEtaOmegaPlus(0),
	
	// - Resolution in Pt as function of generated Pt
    f2dHistAsMCResPtXiMinus(0),		
    f2dHistAsMCResPtXiPlus(0),		
    f2dHistAsMCResPtOmegaMinus(0),
    f2dHistAsMCResPtOmegaPlus(0),	
	
	// - Resolution in R(2D) as function of generated R
    f2dHistAsMCResRXiMinus(0),		
    f2dHistAsMCResRXiPlus(0),		
    f2dHistAsMCResROmegaMinus(0),
    f2dHistAsMCResROmegaPlus(0),

        // - Resolution in phi as function of generated Pt
    f2dHistAsMCResPhiXiMinus(0),
    f2dHistAsMCResPhiXiPlus(0),
    f2dHistAsMCResPhiOmegaMinus(0),
    f2dHistAsMCResPhiOmegaPlus(0),

       //  - Correlation between proton (antiproton) daughter MC pt and Xi/Omega MC pt (to apply Geat/Fluka correction)
    f2dHistAsMCptProtonMCptXiMinus(0),
    f2dHistAsMCptAntiprotonMCptXiPlus(0),
    f2dHistAsMCptProtonMCptOmegaMinus(0),
    f2dHistAsMCptAntiprotonMCptOmegaPlus(0),
    
    fHistV0toXiCosineOfPointingAngle(0),
    fHistV0CosineOfPointingAnglevsPtXi(0),
    fHistV0CosineOfPointingAnglevsPtOmega(0), 
                           
    fCFContCascadePIDAsXiMinus(0),
    fCFContCascadePIDAsXiPlus(0),
    fCFContCascadePIDAsOmegaMinus(0),
    fCFContCascadePIDAsOmegaPlus(0),

    fCFContAsCascadeCuts(0),

    fV0Ampl             (0),
    
    fHistEtaBachXiM  (0),
    fHistEtaPosXiM   (0),
    fHistEtaNegXiM   (0)

{
// Dummy constructor
        for(Int_t iV0selIdx   = 0; iV0selIdx   < 7; iV0selIdx++   ) { fV0VertexerSels          [iV0selIdx   ] = -1.; }
        for(Int_t iCascSelIdx = 0; iCascSelIdx < 8; iCascSelIdx++ ) { fCascadeVertexerSels        [iCascSelIdx ] = -1.; }
}
     
       
     
     
//_____Non-default Constructor________________________________________________________________
AliAnalysisTaskCheckPerformanceCascadePbPb::AliAnalysisTaskCheckPerformanceCascadePbPb(const char *name) 
  : AliAnalysisTaskSE(name),
    fAnalysisType("ESD"), fESDtrackCuts(0), /*fPaveTextBookKeeping(0),*/

    fPIDResponse                   (0),
    fkRerunV0CascVertexers         (0),
    fkQualityCutZprimVtxPos        (kTRUE),
    fkRejectEventPileUp            (kTRUE),
    fkQualityCutNoTPConlyPrimVtx   (kTRUE),
    fkQualityCutTPCrefit           (kTRUE),
    fkQualityCutnTPCcls            (kTRUE),
    fMinnTPCcls                    (0),
    fkExtraSelections              (0),
    fCentrLowLim(0), fCentrUpLim(0), fCentrEstimator(0), fkUseCleaning(0),
    fVtxRange                      (0),
    fApplyAccCut                   (0),
    fMinPtCutOnDaughterTracks      (0),
    fEtaCutOnDaughterTracks      (0),
       
    	// - Cascade part initialisation
    fListHistCascade(0),

    // Events in centraity bins
    fHistEvtsInCentralityBinsvsNtracks(0),
    fHistBestVtxX(0),
    fHistBestVtxY(0),
    fHistBestVtxZ(0),

    // Cascade multiplicity histos

    fHistnXiPlusPerEvTot(0),
    fHistnXiMinusPerEvTot(0),
    fHistnOmegaPlusPerEvTot(0),
    fHistnOmegaMinusPerEvTot(0),

    fHistnXiPlusPerEv(0),
    fHistnXiMinusPerEv(0),
    fHistnOmegaPlusPerEv(0),
    fHistnOmegaMinusPerEv(0),

    fHistnAssoXiMinus(0),
    fHistnAssoXiPlus(0),
    fHistnAssoOmegaMinus(0),
    fHistnAssoOmegaPlus(0),

    fHistMCTrackMultiplicity(0), 
       // - Resolution of the multiplicity estimator
    f2dHistRecoMultVsMCMult(0),

    fHistEtaGenProton(0),
    fHistEtaGenAntiProton(0),
   
// Xi-
   fHistEtaGenCascXiMinus(0),
   f3dHistGenPtVsGenYvsCentXiMinusNat(0),
   f3dHistGenPtVsGenYvsNtracksXiMinusNat(0),
   f3dHistGenPtVsGenYvsCentXiMinusInj(0),
   f3dHistGenPtVsGenYvsNtracksXiMinusInj(0),
   f3dHistGenPtVsGenctauvsCentXiMinusNat(0),
   f3dHistGenPtVsGenctauvsCentXiMinusInj(0),
   
    fHistThetaGenCascXiMinusNat(0),
    fHistThetaGenCascXiMinusInj(0),
    f2dHistGenPtVsGenYFdblXiMinus(0),
    
    fHistThetaLambdaXiMinus(0), 
    fHistThetaBachXiMinus(0),
    
    fHistThetaMesDghterXiMinus(0), 
    fHistThetaBarDghterXiMinus(0),
    
    fHistPtBachXiMinus(0),
    fHistPtMesDghterXiMinus(0),
    fHistPtBarDghterXiMinus(0),
   
    fHistPtRecBachXiMinus(0),
    fHistPtRecMesDghterXiMinus(0),
    fHistPtRecBarDghterXiMinus(0),
 
   // Xi+
   fHistEtaGenCascXiPlus(0),
  f3dHistGenPtVsGenYvsCentXiPlusNat(0),
  f3dHistGenPtVsGenYvsNtracksXiPlusNat(0),
  f3dHistGenPtVsGenYvsCentXiPlusInj(0),
  f3dHistGenPtVsGenYvsNtracksXiPlusInj(0),
  f3dHistGenPtVsGenctauvsCentXiPlusNat(0),
  f3dHistGenPtVsGenctauvsCentXiPlusInj(0),

    fHistThetaGenCascXiPlusNat(0), 
    fHistThetaGenCascXiPlusInj(0),
    f2dHistGenPtVsGenYFdblXiPlus(0),
    
    fHistThetaLambdaXiPlus(0), 
    fHistThetaBachXiPlus(0),
    
    fHistThetaMesDghterXiPlus(0), 
    fHistThetaBarDghterXiPlus(0),
    
    fHistPtBachXiPlus(0),
    fHistPtMesDghterXiPlus(0),
    fHistPtBarDghterXiPlus(0),
   
   // Omega-
   fHistEtaGenCascOmegaMinus(0),
   f3dHistGenPtVsGenYvsCentOmegaMinusNat(0),
   f3dHistGenPtVsGenYvsNtracksOmegaMinusNat(0),
   f3dHistGenPtVsGenYvsCentOmegaMinusInj(0),
   f3dHistGenPtVsGenYvsNtracksOmegaMinusInj(0),
   f3dHistGenPtVsGenctauvsCentOmegaMinusNat(0),
   f3dHistGenPtVsGenctauvsCentOmegaMinusInj(0),
 
    fHistThetaGenCascOmegaMinusNat(0),
    fHistThetaGenCascOmegaMinusInj(0),
    f2dHistGenPtVsGenYFdblOmegaMinus(0),
    
    fHistThetaLambdaOmegaMinus(0), 
    fHistThetaBachOmegaMinus(0),
    
    fHistThetaMesDghterOmegaMinus(0), 
    fHistThetaBarDghterOmegaMinus(0),
    
    fHistPtBachOmegaMinus(0),
    fHistPtMesDghterOmegaMinus(0),
    fHistPtBarDghterOmegaMinus(0),
   
   // Omega+
   fHistEtaGenCascOmegaPlus(0),
   f3dHistGenPtVsGenYvsCentOmegaPlusNat(0),
   f3dHistGenPtVsGenYvsNtracksOmegaPlusNat(0),
   f3dHistGenPtVsGenYvsCentOmegaPlusInj(0),
   f3dHistGenPtVsGenYvsNtracksOmegaPlusInj(0),
   f3dHistGenPtVsGenctauvsCentOmegaPlusNat(0),
   f3dHistGenPtVsGenctauvsCentOmegaPlusInj(0),
   
    fHistThetaGenCascOmegaPlusNat(0),
    fHistThetaGenCascOmegaPlusInj(0),
    f2dHistGenPtVsGenYFdblOmegaPlus(0),
    
    fHistThetaLambdaOmegaPlus(0), 
    fHistThetaBachOmegaPlus(0),
    
    fHistThetaMesDghterOmegaPlus(0), 
    fHistThetaBarDghterOmegaPlus(0),
    
    fHistPtBachOmegaPlus(0),
    fHistPtMesDghterOmegaPlus(0),
    fHistPtBarDghterOmegaPlus(0),

// Part 2 - Association to MC
	
    fHistMassXiMinus(0),
    fHistMassXiPlus(0),
    fHistMassOmegaMinus(0),
    fHistMassOmegaPlus(0),
    
	// - Effective mass histos with combined PID
    fHistMassWithCombPIDXiMinus(0), fHistMassWithCombPIDXiPlus(0),
    fHistMassWithCombPIDOmegaMinus(0), fHistMassWithCombPIDOmegaPlus(0),

	// - PID Probability versus MC Pt(bachelor track)
    f2dHistPIDprobaKaonVsMCPtBach(0), f2dHistPIDprobaPionVsMCPtBach(0),
    
    	// - Effective mass histos with perfect MC PID on the bachelor
    fHistMassWithMcPIDXiMinus(0), fHistMassWithMcPIDXiPlus(0),
    fHistMassWithMcPIDOmegaMinus(0), fHistMassWithMcPIDOmegaPlus(0),
	
    
	// - Effective mass histos for the cascade candidates associated with MC
    fHistAsMCMassXiMinus(0),		
    fHistAsMCMassXiPlus(0),		
    fHistAsMCMassOmegaMinus(0),
    fHistAsMCMassOmegaPlus(0),
    
	// - Generated Pt Vs generated y, for the cascade candidates associated with MC + Info Comb. PID
    f2dHistAsMCandCombPIDGenPtVsGenYXiMinus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYXiPlus(0),
    f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus(0), 
    f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus(0),
    
    	// - Generated Pt Vs generated y, for the cascade candidates associated with MC
    f2dHistAsMCGenPtVsGenYXiMinus(0),
    f2dHistAsMCGenPtVsGenYXiPlus(0),
    f2dHistAsMCGenPtVsGenYOmegaMinus(0),
    f2dHistAsMCGenPtVsGenYOmegaPlus(0),
    
    	// - Generated Eta of the the cascade candidates associated with MC
    fHistAsMCGenEtaXiMinus(0),
    fHistAsMCGenEtaXiPlus(0),
    fHistAsMCGenEtaOmegaMinus(0),
    fHistAsMCGenEtaOmegaPlus(0),
	
	// - Resolution in Pt as function of generated Pt
    f2dHistAsMCResPtXiMinus(0),		
    f2dHistAsMCResPtXiPlus(0),		
    f2dHistAsMCResPtOmegaMinus(0),
    f2dHistAsMCResPtOmegaPlus(0),	
	
	// - Resolution in R(2D) as function of generated R
    f2dHistAsMCResRXiMinus(0),		
    f2dHistAsMCResRXiPlus(0),		
    f2dHistAsMCResROmegaMinus(0),
    f2dHistAsMCResROmegaPlus(0),

       // - Resolution in phi as function of generated Pt
    f2dHistAsMCResPhiXiMinus(0),
    f2dHistAsMCResPhiXiPlus(0),
    f2dHistAsMCResPhiOmegaMinus(0),
    f2dHistAsMCResPhiOmegaPlus(0),

       //  - Correlation between proton (antiproton) daughter MC pt and Xi/Omega MC pt (to apply Geat/Fluka correction)
    f2dHistAsMCptProtonMCptXiMinus(0),
    f2dHistAsMCptAntiprotonMCptXiPlus(0),
    f2dHistAsMCptProtonMCptOmegaMinus(0),
    f2dHistAsMCptAntiprotonMCptOmegaPlus(0),

    fHistV0toXiCosineOfPointingAngle(0),
    fHistV0CosineOfPointingAnglevsPtXi(0),
    fHistV0CosineOfPointingAnglevsPtOmega(0),

    fCFContCascadePIDAsXiMinus(0),
    fCFContCascadePIDAsXiPlus(0),
    fCFContCascadePIDAsOmegaMinus(0),
    fCFContCascadePIDAsOmegaPlus(0),

    fCFContAsCascadeCuts(0),

    fV0Ampl(0),

    fHistEtaBachXiM  (0),
    fHistEtaPosXiM   (0),
    fHistEtaNegXiM   (0)


{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // Output slot #1 writes into a TList container (cascade)
        
        // - default PbPb cuts
        //fV0VertexerSels[0] =  33.  ;  // max allowed chi2
        //fV0VertexerSels[1] =   0.1; // min allowed impact parameter for the 1st daughter 
        //fV0VertexerSels[2] =   0.1; // min allowed impact parameter for the 2nd daughter 
        //fV0VertexerSels[3] =   1.0 ;  // max allowed DCA between the daughter tracks       
        //fV0VertexerSels[4] =   0.998 ;  // min allowed cosine of V0's pointing angle         
        //fV0VertexerSels[5] =   0.9;  // min radius of the fiducial volume                 
        //fV0VertexerSels[6] = 100.  ;  // max radius of the fiducial volume                 
        //fCascVertexerSels[0] =  33.   ;  // max allowed chi2 
        //fCascVertexerSels[1] =   0.05;  // min allowed V0 impact parameter                    
        //fCascVertexerSels[2] =   0.008;  // "window" around the Lambda mass                    
        //fCascVertexerSels[3] =   0.03;  // min allowed bachelor's impact parameter            
        //fCascVertexerSels[4] =   0.3  ;  // max allowed DCA between the V0 and the bachelor    
        //fCascVertexerSels[5] =   0.999;  // min allowed cosine of the cascade pointing angle   
        //fCascVertexerSels[6] =   0.9  ;  // min radius of the fiducial volume                  
        //fCascVertexerSels[7] = 100.   ;  // max radius of the fiducial volume                  
        // - default pp values
        fV0VertexerSels[0] =  33.  ;        // max allowed chi2
        fV0VertexerSels[1] =   0.05;        // min allowed impact parameter for the 1st daughter 
        fV0VertexerSels[2] =   0.05;        // min allowed impact parameter for the 2nd daughter 
        fV0VertexerSels[3] =   1.5 ;        // max allowed DCA between the daughter tracks       
        fV0VertexerSels[4] =   0.9 ;        // min allowed cosine of V0's pointing angle         
        fV0VertexerSels[5] =   0.2 ;        // min radius of the fiducial volume                 
        fV0VertexerSels[6] = 100.  ;        // max radius of the fiducial volume                 
        fCascadeVertexerSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
        fCascadeVertexerSels[1] =   0.01 ;  // min allowed V0 impact parameter                    
        fCascadeVertexerSels[2] =   0.008;  // "window" around the Lambda mass                    
        fCascadeVertexerSels[3] =   0.01 ;  // min allowed bachelor's impact parameter 
        fCascadeVertexerSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    
        fCascadeVertexerSels[5] =   0.98 ;  // min allowed cosine of the cascade pointing angle   
        fCascadeVertexerSels[6] =   0.9  ;  // min radius of the fiducial volume                  
        fCascadeVertexerSels[7] = 100.   ;  // max radius of the fiducial volume 
        
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliCFContainer::Class());
  DefineOutput(3, AliCFContainer::Class());
  DefineOutput(4, AliCFContainer::Class());
  DefineOutput(5, AliCFContainer::Class());
  DefineOutput(6, AliCFContainer::Class()); 

}


AliAnalysisTaskCheckPerformanceCascadePbPb::~AliAnalysisTaskCheckPerformanceCascadePbPb()
{
  //
  // Destructor
  //

  // For all TH1, 2, 3 HnSparse and CFContainer are in the fListCascade TList.
  // They will be deleted when fListCascade is deleted by the TSelector dtor
  // Because of TList::SetOwner()

  if (fListHistCascade && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())      { delete fListHistCascade;     fListHistCascade = 0x0;  }  
  if (fCFContCascadePIDAsXiMinus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fCFContCascadePIDAsXiMinus;     fCFContCascadePIDAsXiMinus = 0x0;  }
  if (fCFContCascadePIDAsXiPlus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fCFContCascadePIDAsXiPlus;     fCFContCascadePIDAsXiPlus = 0x0;  }
  if (fCFContCascadePIDAsOmegaMinus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fCFContCascadePIDAsOmegaMinus;     fCFContCascadePIDAsOmegaMinus = 0x0;  }
  if (fCFContCascadePIDAsOmegaPlus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fCFContCascadePIDAsOmegaPlus;     fCFContCascadePIDAsOmegaPlus = 0x0;  }
  if (fCFContAsCascadeCuts && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fCFContAsCascadeCuts;     fCFContAsCascadeCuts = 0x0;  }

  if (fESDtrackCuts)         { delete fESDtrackCuts;        fESDtrackCuts = 0x0; }
  /*if (fPaveTextBookKeeping)  { delete fPaveTextBookKeeping; fPaveTextBookKeeping = 0x0; } // fPaveTextBookKeeping is not stored into the TList*/
}


//________________________________________________________________________
void AliAnalysisTaskCheckPerformanceCascadePbPb::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

   // Option for AliLog
	AliLog::SetGlobalLogLevel(AliLog::kError); 
   	// to suppress the extensive info prompted by a run with MC			

   // Definition of the output datamembers	
   fListHistCascade = new TList();
   fListHistCascade->SetOwner(); // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

   // New PID object
   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
   fPIDResponse = inputHandler->GetPIDResponse();
        
   // Only used to get the number of primary reconstructed tracks
   if (! fESDtrackCuts ){
     fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE); // Std definition of primary (see kTRUE argument) tracks for 2010
//     fESDtrackCuts->SetEtaRange(-0.8,+0.8);
//     fESDtrackCuts->SetPtRange(0.15, 1e10);
     Printf("CheckCascade - ESDtrackCuts set up to 2010 std ITS-TPC cuts...");
   }


/*
if( !fPaveTextBookKeeping){
        fPaveTextBookKeeping = new TPaveText(0.1, 0.1, 0.9, 0.9,"NDC");
        fPaveTextBookKeeping->SetName("fPaveTextBookKeeping");
        fPaveTextBookKeeping->SetBorderSize(0);
        fPaveTextBookKeeping->SetTextAlign(12);
        fPaveTextBookKeeping->SetFillColor(kWhite);
        fPaveTextBookKeeping->SetTextFont(42);        // regular Arial or Helvetica,
        fPaveTextBookKeeping->SetTextColor(kGray+3);
        
        
        fPaveTextBookKeeping->AddText( "Task CHECK PERFORMANCE CASCADE analysis" );
        fPaveTextBookKeeping->AddText("- - - - - - - - - - - ");
        fPaveTextBookKeeping->AddText( Form("AnalysisType : %s ", fAnalysisType.Data() ));
        fPaveTextBookKeeping->AddText("Colliding system : A-A collisions ");

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
        if(fkQualityCutnTPCcls)          fPaveTextBookKeeping->AddText("E. Quality Cut(n TPC clusters)   = ON  ");
        else                             fPaveTextBookKeeping->AddText("E. Quality Cut(n TPC clusters)   = Off ");
        if(fkExtraSelections)            fPaveTextBookKeeping->AddText("F. Extra Analysis Selections         = ON  ");
        else                             fPaveTextBookKeeping->AddText("F. Extra Analysis Selections         = Off ");

        fPaveTextBookKeeping->AddText("- - - - - - - - - - - ");

        fListHistCascade->Add(fPaveTextBookKeeping);
}       
*/
 
  // - General
  Double_t ptBinLimits[101];
  for (Int_t iptbin = 0; iptbin<101; ++iptbin) ptBinLimits[iptbin]=iptbin*0.1;
  Double_t yBinLimits[111];
  for (Int_t iybin = 0; iybin<111; ++iybin) yBinLimits[iybin]=-1.1+iybin*0.02;

  // Events in centrality bins
  Double_t centBinLimits[12] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};

  Double_t ctauBinLimits[112];
  for (Int_t ict = 0; ict<112; ++ict) ctauBinLimits[ict] = (Double_t) (ict-1.); 

  fHistEvtsInCentralityBinsvsNtracks = new TH2F("fHistEvtsInCentralityBinsvsNtracks","",11,centBinLimits,100,0.,4000.);
  fListHistCascade->Add(fHistEvtsInCentralityBinsvsNtracks);

  fHistBestVtxX = new TH1F("fHistBestVtxX","",100,-0.1,0.1);
  fListHistCascade->Add(fHistBestVtxX);
  fHistBestVtxY = new TH1F("fHistBestVtxY","",100,-1,1);
  fListHistCascade->Add(fHistBestVtxY);
  fHistBestVtxZ = new TH1F("fHistBestVtxZ","",100,-20.,20.);
  fListHistCascade->Add(fHistBestVtxZ);

  // Cascade multiplicity distributions
  
  fHistnXiPlusPerEvTot= new TH1F("fHistnXiPlusPerEvTot", "", 100, 0, 100);
  fListHistCascade->Add(fHistnXiPlusPerEvTot);
  fHistnXiMinusPerEvTot= new TH1F("fHistnXiMinusPerEvTot", "", 100, 0, 100);
  fListHistCascade->Add(fHistnXiMinusPerEvTot);
  fHistnOmegaPlusPerEvTot = new TH1F("fHistnOmegaPlusPerEvTot", "", 50, 0, 50);
  fListHistCascade->Add(fHistnOmegaPlusPerEvTot);
  fHistnOmegaMinusPerEvTot= new TH1F("fHistnOmegaMinusPerEvTot", "", 50, 0, 50);
  fListHistCascade->Add(fHistnOmegaMinusPerEvTot);
     
  fHistnXiPlusPerEv= new TH1F("fHistnXiPlusPerEv", "", 100, 0, 100);
  fListHistCascade->Add(fHistnXiPlusPerEv);
  fHistnXiMinusPerEv= new TH1F("fHistnXiMinusPerEv", "", 100, 0, 100);
  fListHistCascade->Add(fHistnXiMinusPerEv);
  fHistnOmegaPlusPerEv= new TH1F("fHistnOmegaPlusPerEv", "", 50, 0, 50);
  fListHistCascade->Add(fHistnOmegaPlusPerEv);
  fHistnOmegaMinusPerEv= new TH1F("fHistnOmegaMinusPerEv", "", 50, 0, 50);
  fListHistCascade->Add(fHistnOmegaMinusPerEv);

  fHistnAssoXiMinus= new TH1F("fHistnAssoXiMinus", "", 100, 0, 100);
  fListHistCascade->Add(fHistnAssoXiMinus);
  fHistnAssoXiPlus= new TH1F("fHistnAssoXiPlus", "", 100, 0, 100);
  fListHistCascade->Add(fHistnAssoXiPlus); 
  fHistnAssoOmegaMinus= new TH1F("fHistnAssoOmegaMinus", "", 50, 0, 50);
  fListHistCascade->Add(fHistnAssoOmegaMinus);
  fHistnAssoOmegaPlus= new TH1F("fHistnAssoOmegaPlus", "", 50, 0, 50);
  fListHistCascade->Add(fHistnAssoOmegaPlus);
  
  if (!fHistMCTrackMultiplicity) {
    fHistMCTrackMultiplicity = new TH1F("fHistMCTrackMultiplicity", "Multiplicity distribution;Number of tracks;Events", 200, 0, 200000.); 
    fListHistCascade->Add(fHistMCTrackMultiplicity);
  }
  
    // - Resolution of the multiplicity estimator
  if (! f2dHistRecoMultVsMCMult){
    f2dHistRecoMultVsMCMult = new TH2F("f2dHistRecoMultVsMCMult", "Resolution of the multiplicity estimator (prim. tracks); Reco Multiplicity (prim. tracks); MC multiplicity (gen. part.)", 200, 0., 150000., 200, 0., 200000.);
    fListHistCascade->Add(f2dHistRecoMultVsMCMult);
  }
  
  if (!fHistEtaGenProton) {
    fHistEtaGenProton = new TH1F("fHistEtaGenProton", "#eta of any gen. p^{+};#eta;Number of prim. protons", 200, -10, 10);
    fListHistCascade->Add(fHistEtaGenProton);
  }
  
  if (!fHistEtaGenAntiProton) {
     fHistEtaGenAntiProton = new TH1F("fHistEtaGenAntiProton", "#eta of any gen. #bar{p}^{-};#eta;Number of prim. #bar{p}", 200, -10, 10);
     fListHistCascade->Add(fHistEtaGenAntiProton);
  }
 
 
  //--------
  // I - Xi- 
  // - Pseudo-Rapidity distribution
  if (!fHistEtaGenCascXiMinus) {
     fHistEtaGenCascXiMinus = new TH1F("fHistEtaGenCascXiMinus", "#eta of any gen. #Xi^{-};#eta;Number of Casc", 200, -10, 10);
     fListHistCascade->Add(fHistEtaGenCascXiMinus);
   }

  if (!f3dHistGenPtVsGenYvsCentXiMinusNat) {
     f3dHistGenPtVsGenYvsCentXiMinusNat = new TH3D("f3dHistGenPtVsGenYvsCentXiMinusNat", "MC P_{t} Vs MC Y of Gen #Xi^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 100, ptBinLimits, 110, yBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsCentXiMinusNat);
   }
  if (!f3dHistGenPtVsGenYvsNtracksXiMinusNat) {
    f3dHistGenPtVsGenYvsNtracksXiMinusNat = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinusNat", "MC P_{t} Vs MC Y of Gen #Xi^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 100, 0., 4000.);
    fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinusNat);
  }
  if (!f3dHistGenPtVsGenYvsCentXiMinusInj) {
     f3dHistGenPtVsGenYvsCentXiMinusInj = new TH3D("f3dHistGenPtVsGenYvsCentXiMinusInj", "MC P_{t} Vs MC Y of Gen #Xi^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 100, ptBinLimits, 110, yBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsCentXiMinusInj);
  }
  if (!f3dHistGenPtVsGenYvsNtracksXiMinusInj) {
    f3dHistGenPtVsGenYvsNtracksXiMinusInj = new TH3D("f3dHistGenPtVsGenYvsNtracksXiMinusInj", "MC P_{t} Vs MC Y of Gen #Xi^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 100, 0., 4000.);
    fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiMinusInj);
  }
  if (!f3dHistGenPtVsGenctauvsCentXiMinusNat) {
     f3dHistGenPtVsGenctauvsCentXiMinusNat = new TH3D("f3dHistGenPtVsGenctauvsCentXiMinusNat", "MC P_{t} Vs MC ctau Vs Centrality of Gen #Xi^{-} ", 100, ptBinLimits, 111, ctauBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenctauvsCentXiMinusNat);
  }

  if (!f3dHistGenPtVsGenctauvsCentXiMinusInj) {
     f3dHistGenPtVsGenctauvsCentXiMinusInj = new TH3D("f3dHistGenPtVsGenctauvsCentXiMinusInj", "MC P_{t} Vs MC ctau Vs Centrality of Gen #Xi^{-} ", 100, ptBinLimits, 111, ctauBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenctauvsCentXiMinusInj);
  }
  

  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascXiMinusNat) {
     fHistThetaGenCascXiMinusNat = new TH1F("fHistThetaGenCascXiMinusNat", "#theta of gen. #Xi^{-};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascXiMinusNat);
  }

  if (!fHistThetaGenCascXiMinusInj) {
     fHistThetaGenCascXiMinusInj = new TH1F("fHistThetaGenCascXiMinusInj", "#theta of injected. #Xi^{-};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascXiMinusInj);
  }


  if (!f2dHistGenPtVsGenYFdblXiMinus) {
     f2dHistGenPtVsGenYFdblXiMinus = new TH2D("f2dHistGenPtVsGenYFdblXiMinus", "MC P_{t} Vs MC Y of findable Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistGenPtVsGenYFdblXiMinus);
  }
  
		// - Theta distribution the daughters (control plots)
  
  if (!fHistThetaLambdaXiMinus) {
     fHistThetaLambdaXiMinus = new TH1F("fHistThetaLambdaXiMinus", "#theta of gen. #Lambda (Xi dghter);#theta_{#Lambda};Number of #Lambda^0", 200, -10, 190);
     fListHistCascade->Add(fHistThetaLambdaXiMinus);
  }

  if (!fHistThetaBachXiMinus) {
     fHistThetaBachXiMinus = new TH1F("fHistThetaBachXiMinus", "#theta of gen. Bach.;#theta_{Bach};Number of Bach.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaBachXiMinus);
  }
  
  if (!fHistThetaMesDghterXiMinus) {
     fHistThetaMesDghterXiMinus = new TH1F("fHistThetaMesDghterXiMinus", "#theta of gen. Meson #Lambda dghter;#theta_{MesDght};Number of Mes.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaMesDghterXiMinus);
  }
  
  if (!fHistThetaBarDghterXiMinus) {
     fHistThetaBarDghterXiMinus = new TH1F("fHistThetaBarDghterXiMinus", "#theta of gen. Baryon #Lambda dghter;#theta_{BarDght};Number of Bar.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaBarDghterXiMinus);
  }
 
  		// - Pt distribution (control plots)
    
  if (!fHistPtBachXiMinus) {
     fHistPtBachXiMinus = new TH1F("fHistPtBachXiMinus", "p_{t} of gen. Bach.;pt_{Bach};Number of Bach.", 200, 0, 10);
     fListHistCascade->Add(fHistPtBachXiMinus);
  }
  
  if (!fHistPtMesDghterXiMinus) {
     fHistPtMesDghterXiMinus = new TH1F("fHistPtMesDghterXiMinus", "p_{t} of gen. Meson #Lambda dghter;pt_{MesDght};Number of Mes.", 200, 0, 10);
     fListHistCascade->Add(fHistPtMesDghterXiMinus);
  }
    
  if (!fHistPtBarDghterXiMinus) {
     fHistPtBarDghterXiMinus = new TH1F("fHistPtBarDghterXiMinus", "p_{t} of gen. Baryon #Lambda dghter;pt_{BarDght};Number of Bar.", 200, 0, 10);
     fListHistCascade->Add(fHistPtBarDghterXiMinus);
  }
  
  if (!fHistPtRecBachXiMinus) {
     fHistPtRecBachXiMinus = new TH1F("fHistPtRecBachXiMinus", "p_{t} of Bach.;pt_{Bach};Number of Bach.", 200, 0, 10);
     fListHistCascade->Add(fHistPtRecBachXiMinus);
  }

  if (!fHistPtRecMesDghterXiMinus) {
     fHistPtRecMesDghterXiMinus = new TH1F("fHistPtRecMesDghterXiMinus", "p_{t} of Meson #Lambda dghter;pt_{MesDght};Number of Mes.", 200, 0, 10);
     fListHistCascade->Add(fHistPtRecMesDghterXiMinus);
  }

  if (!fHistPtRecBarDghterXiMinus) {
     fHistPtRecBarDghterXiMinus = new TH1F("fHistPtRecBarDghterXiMinus", "p_{t} of Baryon #Lambda dghter;pt_{BarDght};Number of Bar.", 200, 0, 10);
     fListHistCascade->Add(fHistPtRecBarDghterXiMinus);
  }
 


 
  
  //--------
  // II - Xi+ 
  // - Pseudo-Rapidity distribution
  if (!fHistEtaGenCascXiPlus) {
     fHistEtaGenCascXiPlus = new TH1F("fHistEtaGenCascXiPlus", "#eta of any gen. #Xi^{+};#eta;Number of Casc", 200, -10, 10);
     fListHistCascade->Add(fHistEtaGenCascXiPlus);
  }
  if (!f3dHistGenPtVsGenYvsCentXiPlusNat) {
     f3dHistGenPtVsGenYvsCentXiPlusNat = new TH3D("f3dHistGenPtVsGenYvsCentXiPlusNat", "MC P_{t} Vs MC Y of Gen #Xi^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 100, ptBinLimits, 110, yBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsCentXiPlusNat);
  }
  if (!f3dHistGenPtVsGenYvsNtracksXiPlusNat) {
    f3dHistGenPtVsGenYvsNtracksXiPlusNat = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlusNat", "MC P_{t} Vs MC Y of Gen #Xi^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 100, 0., 4000.);
    fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlusNat);
  }
  if (!f3dHistGenPtVsGenYvsCentXiPlusInj) {
     f3dHistGenPtVsGenYvsCentXiPlusInj = new TH3D("f3dHistGenPtVsGenYvsCentXiPlusInj", "MC P_{t} Vs MC Y of Gen #Xi^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 100, ptBinLimits, 110, yBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsCentXiPlusInj);
  }
  if (!f3dHistGenPtVsGenYvsNtracksXiPlusInj) {
    f3dHistGenPtVsGenYvsNtracksXiPlusInj = new TH3D("f3dHistGenPtVsGenYvsNtracksXiPlusInj", "MC P_{t} Vs MC Y of Gen #Xi^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 100, 0., 4000.);
    fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksXiPlusInj);
  }
  if (!f3dHistGenPtVsGenctauvsCentXiPlusNat) {
     f3dHistGenPtVsGenctauvsCentXiPlusNat = new TH3D("f3dHistGenPtVsGenctauvsCentXiPlusNat", "MC P_{t} Vs MC ctau Vs Centrality of Gen #Xi^{+} ", 100, ptBinLimits, 111, ctauBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenctauvsCentXiPlusNat);
  }
  if (!f3dHistGenPtVsGenctauvsCentXiPlusInj) {
     f3dHistGenPtVsGenctauvsCentXiPlusInj = new TH3D("f3dHistGenPtVsGenctauvsCentXiPlusInj", "MC P_{t} Vs MC ctau Vs Centrality of Gen #Xi^{+} ", 100, ptBinLimits, 111, ctauBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenctauvsCentXiPlusInj);
  }

  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascXiPlusNat) {
     fHistThetaGenCascXiPlusNat = new TH1F("fHistThetaGenCascXiPlusNat", "#theta of gen. #Xi^{+};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascXiPlusNat);
  }
  if (!fHistThetaGenCascXiPlusInj) {
     fHistThetaGenCascXiPlusInj = new TH1F("fHistThetaGenCascXiPlusInj", "#theta of inj. #Xi^{+};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascXiPlusInj);
  }

 
  if (!f2dHistGenPtVsGenYFdblXiPlus) {
     f2dHistGenPtVsGenYFdblXiPlus = new TH2D("f2dHistGenPtVsGenYFdblXiPlus", "MC P_{t} Vs MC Y of findable Gen #Xi^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistGenPtVsGenYFdblXiPlus);
  }
  
  		// - Theta distribution the daughters (control plots)
  
  if (!fHistThetaLambdaXiPlus) {
     fHistThetaLambdaXiPlus = new TH1F("fHistThetaLambdaXiPlus", "#theta of gen. #Lambda (Xi dghter);#theta_{#Lambda};Number of #Lambda", 200, -10, 190);
     fListHistCascade->Add(fHistThetaLambdaXiPlus);
  }

  if (!fHistThetaBachXiPlus) {
     fHistThetaBachXiPlus = new TH1F("fHistThetaBachXiPlus", "#theta of gen. Bach.;#theta_{Bach};Number of Bach.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaBachXiPlus);
  }
  
  if (!fHistThetaMesDghterXiPlus) {
     fHistThetaMesDghterXiPlus = new TH1F("fHistThetaMesDghterXiPlus", "#theta of gen. Meson #Lambda dghter;#theta_{MesDght};Number of Mes.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaMesDghterXiPlus);
  }
  
  if (!fHistThetaBarDghterXiPlus) {
     fHistThetaBarDghterXiPlus = new TH1F("fHistThetaBarDghterXiPlus", "#theta of gen. Baryon #Lambda dghter;#theta_{BarDght};Number of Bar.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaBarDghterXiPlus);
  }
 
  		// - Pt distribution (control plots)
    
  if (!fHistPtBachXiPlus) {
     fHistPtBachXiPlus = new TH1F("fHistPtBachXiPlus", "p_{t} of gen. Bach.;pt_{Bach};Number of Bach.", 200, 0, 10);
     fListHistCascade->Add(fHistPtBachXiPlus);
  }
  
  if (!fHistPtMesDghterXiPlus) {
     fHistPtMesDghterXiPlus = new TH1F("fHistPtMesDghterXiPlus", "p_{t} of gen. Meson #Lambda dghter);pt_{MesDght};Number of Mes.", 200, 0, 10);
     fListHistCascade->Add(fHistPtMesDghterXiPlus);
  }
    
  if (!fHistPtBarDghterXiPlus) {
     fHistPtBarDghterXiPlus = new TH1F("fHistPtBarDghterXiPlus", "p_{t} of gen. Baryon #Lambda dghter);pt_{BarDght};Number of Bar.", 200, 0, 10);
     fListHistCascade->Add(fHistPtBarDghterXiPlus);
  }
  
  
  //---------
  // III - Omega- 
  // - Pseudo-Rapidity distribution
  if (!fHistEtaGenCascOmegaMinus) {
     fHistEtaGenCascOmegaMinus = new TH1F("fHistEtaGenCascOmegaMinus", "#eta of any gen. #Omega^{-};#eta;Number of Casc", 200, -10, 10);
     fListHistCascade->Add(fHistEtaGenCascOmegaMinus);
  }
 if (!f3dHistGenPtVsGenYvsCentOmegaMinusNat) {
     f3dHistGenPtVsGenYvsCentOmegaMinusNat = new TH3D("f3dHistGenPtVsGenYvsCentOmegaMinusNat", "MC P_{t} Vs MC Y of Gen #Omega^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 100, ptBinLimits, 110, yBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsCentOmegaMinusNat);
  }
  if (!f3dHistGenPtVsGenYvsNtracksOmegaMinusNat) {
     f3dHistGenPtVsGenYvsNtracksOmegaMinusNat = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinusNat", "MC P_{t} Vs MC Y of Gen #Omega^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 100, 0., 4000.);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinusNat);
  }
  if (!f3dHistGenPtVsGenYvsCentOmegaMinusInj) {
     f3dHistGenPtVsGenYvsCentOmegaMinusInj = new TH3D("f3dHistGenPtVsGenYvsCentOmegaMinusInj", "MC P_{t} Vs MC Y of Gen #Omega^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 100, ptBinLimits, 110, yBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsCentOmegaMinusInj);
  }
  if (!f3dHistGenPtVsGenYvsNtracksOmegaMinusInj) {
     f3dHistGenPtVsGenYvsNtracksOmegaMinusInj = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaMinusInj", "MC P_{t} Vs MC Y of Gen #Omega^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 100, 0., 4000.);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaMinusInj);
  }
  if (!f3dHistGenPtVsGenctauvsCentOmegaMinusNat) { 
     f3dHistGenPtVsGenctauvsCentOmegaMinusNat = new TH3D("f3dHistGenPtVsGenctauvsCentOmegaMinusNat", "MC P_{t} Vs MC ctau Vs Centrality of Gen #Omega^{-} ", 100, ptBinLimits, 111, ctauBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenctauvsCentOmegaMinusNat);
  }
  if (!f3dHistGenPtVsGenctauvsCentOmegaMinusInj) {
     f3dHistGenPtVsGenctauvsCentOmegaMinusInj = new TH3D("f3dHistGenPtVsGenctauvsCentOmegaMinusInj", "MC P_{t} Vs MC ctau Vs Centrality of Gen #Omega^{-} ", 100, ptBinLimits, 111, ctauBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenctauvsCentOmegaMinusInj);
  }



  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascOmegaMinusNat) {
     fHistThetaGenCascOmegaMinusNat = new TH1F("fHistThetaGenCascOmegaMinusNat", "#theta of gen. #Omega^{-};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascOmegaMinusNat);
  }
  if (!fHistThetaGenCascOmegaMinusInj) {
     fHistThetaGenCascOmegaMinusInj = new TH1F("fHistThetaGenCascOmegaMinusInj", "#theta of inj. #Omega^{-};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascOmegaMinusInj);
  }

 
  if (!f2dHistGenPtVsGenYFdblOmegaMinus) {
     f2dHistGenPtVsGenYFdblOmegaMinus = new TH2D("f2dHistGenPtVsGenYFdblOmegaMinus", "MC P_{t} Vs MC Y of findable Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistGenPtVsGenYFdblOmegaMinus);
  }
  
  		// - Theta distribution the daughters (control plots)
  
  if (!fHistThetaLambdaOmegaMinus) {
     fHistThetaLambdaOmegaMinus = new TH1F("fHistThetaLambdaOmegaMinus", "#theta of gen. #Lambda (Omega dghter);#theta_{#Lambda};Number of #Lambda", 200, -10, 190);
     fListHistCascade->Add(fHistThetaLambdaOmegaMinus);
  }

  if (!fHistThetaBachOmegaMinus) {
     fHistThetaBachOmegaMinus = new TH1F("fHistThetaBachOmegaMinus", "#theta of gen. Bach.;#theta_{Bach};Number of Bach.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaBachOmegaMinus);
  }
  
  if (!fHistThetaMesDghterOmegaMinus) {
     fHistThetaMesDghterOmegaMinus = new TH1F("fHistThetaMesDghterOmegaMinus", "#theta of gen. Meson #Lambda dghter;#theta_{MesDght};Number of Mes.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaMesDghterOmegaMinus);
  }
  
  if (!fHistThetaBarDghterOmegaMinus) {
     fHistThetaBarDghterOmegaMinus = new TH1F("fHistThetaBarDghterOmegaMinus", "#theta of gen. Baryon #Lambda dghter;#theta_{BarDght};Number of Bar.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaBarDghterOmegaMinus);
  }
 
  		// - Pt distribution (control plots)
    
  if (!fHistPtBachOmegaMinus) {
     fHistPtBachOmegaMinus = new TH1F("fHistPtBachOmegaMinus", "p_{t} of gen. Bach.;pt_{Bach};Number of Bach.", 200, 0, 10);
     fListHistCascade->Add(fHistPtBachOmegaMinus);
  }
  
  if (!fHistPtMesDghterOmegaMinus) {
     fHistPtMesDghterOmegaMinus = new TH1F("fHistPtMesDghterOmegaMinus", "p_{t} of gen. Meson #Lambda dghter);pt_{MesDght};Number of Mes.", 200, 0, 10);
     fListHistCascade->Add(fHistPtMesDghterOmegaMinus);
  }
    
  if (!fHistPtBarDghterOmegaMinus) {
     fHistPtBarDghterOmegaMinus = new TH1F("fHistPtBarDghterOmegaMinus", "p_{t} of gen. Baryon #Lambda dghter);pt_{BarDght};Number of Bar.", 200, 0, 10);
     fListHistCascade->Add(fHistPtBarDghterOmegaMinus);
  }
  
  
  //---------
  // IV - Omega+ 
  // - Pseudo-Rapidity distribution
  if (!fHistEtaGenCascOmegaPlus) {
     fHistEtaGenCascOmegaPlus = new TH1F("fHistEtaGenCascOmegaPlus", "#eta of any gen. #Omega^{+};#eta;Number of Casc", 200, -10, 10);
     fListHistCascade->Add(fHistEtaGenCascOmegaPlus);
  }
  if (!f3dHistGenPtVsGenYvsCentOmegaPlusNat) {
     f3dHistGenPtVsGenYvsCentOmegaPlusNat = new TH3D("f3dHistGenPtVsGenYvsCentOmegaPlusNat", "MC P_{t} Vs MC Y of Gen #Omega^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 100, ptBinLimits, 110, yBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsCentOmegaPlusNat);
  }
  if (!f3dHistGenPtVsGenYvsNtracksOmegaPlusNat) {
     f3dHistGenPtVsGenYvsNtracksOmegaPlusNat = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlusNat", "MC P_{t} Vs MC Y of Gen #Omega^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 100, 0., 4000.);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlusNat);
  }
  if (!f3dHistGenPtVsGenYvsCentOmegaPlusInj) {
     f3dHistGenPtVsGenYvsCentOmegaPlusInj = new TH3D("f3dHistGenPtVsGenYvsCentOmegaPlusInj", "MC P_{t} Vs MC Y of Gen #Omega^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 100, ptBinLimits, 110, yBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsCentOmegaPlusInj);
  }
  if (!f3dHistGenPtVsGenYvsNtracksOmegaPlusInj) {
     f3dHistGenPtVsGenYvsNtracksOmegaPlusInj = new TH3D("f3dHistGenPtVsGenYvsNtracksOmegaPlusInj", "MC P_{t} Vs MC Y of Gen #Omega^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 110, -1.1, 1.1, 100, 0., 4000.);
     fListHistCascade->Add(f3dHistGenPtVsGenYvsNtracksOmegaPlusInj);
  }
  if (!f3dHistGenPtVsGenctauvsCentOmegaPlusNat) {
     f3dHistGenPtVsGenctauvsCentOmegaPlusNat = new TH3D("f3dHistGenPtVsGenctauvsCentOmegaPlusNat", "MC P_{t} Vs MC ctau Vs Centrality of Gen #Omega^{+} ", 100, ptBinLimits, 111, ctauBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenctauvsCentOmegaPlusNat);
  }
  if (!f3dHistGenPtVsGenctauvsCentOmegaPlusInj) {
     f3dHistGenPtVsGenctauvsCentOmegaPlusInj = new TH3D("f3dHistGenPtVsGenctauvsCentOmegaPlusInj", "MC P_{t} Vs MC ctau Vs Centrality of Gen #Omega^{+} ", 100, ptBinLimits, 111, ctauBinLimits, 11, centBinLimits);
     fListHistCascade->Add(f3dHistGenPtVsGenctauvsCentOmegaPlusInj);
  }
 
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascOmegaPlusNat) {
     fHistThetaGenCascOmegaPlusNat = new TH1F("fHistThetaGenCascOmegaPlusNat", "#theta of gen. #Omega^{+};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascOmegaPlusNat);
  }
  if (!fHistThetaGenCascOmegaPlusInj) {
     fHistThetaGenCascOmegaPlusInj = new TH1F("fHistThetaGenCascOmegaPlusInj", "#theta of inj. #Omega^{+};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascOmegaPlusInj);
  }

 
  if (!f2dHistGenPtVsGenYFdblOmegaPlus) {
     f2dHistGenPtVsGenYFdblOmegaPlus = new TH2D("f2dHistGenPtVsGenYFdblOmegaPlus", "MC P_{t} Vs MC Y of findable Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 100, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistGenPtVsGenYFdblOmegaPlus);
  }

  
  		// - Theta distribution the daughters (control plots)
  
  if (!fHistThetaLambdaOmegaPlus) {
     fHistThetaLambdaOmegaPlus = new TH1F("fHistThetaLambdaOmegaPlus", "#theta of gen. #Lambda (Omega dghter);#theta_{#Lambda};Number of #Lambda", 200, -10, 190);
     fListHistCascade->Add(fHistThetaLambdaOmegaPlus);
  }

  if (!fHistThetaBachOmegaPlus) {
     fHistThetaBachOmegaPlus = new TH1F("fHistThetaBachOmegaPlus", "#theta of gen. Bach.;#theta_{Bach};Number of Bach.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaBachOmegaPlus);
  }
  
  if (!fHistThetaMesDghterOmegaPlus) {
     fHistThetaMesDghterOmegaPlus = new TH1F("fHistThetaMesDghterOmegaPlus", "#theta of gen. Meson #Lambda dghter;#theta_{MesDght};Number of Mes.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaMesDghterOmegaPlus);
  }
  
  if (!fHistThetaBarDghterOmegaPlus) {
     fHistThetaBarDghterOmegaPlus = new TH1F("fHistThetaBarDghterOmegaPlus", "#theta of gen. Baryon #Lambda dghter;#theta_{BarDght};Number of Bar.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaBarDghterOmegaPlus);
  }
 
  		// - Pt distribution (control plots)
    
  if (!fHistPtBachOmegaPlus) {
     fHistPtBachOmegaPlus = new TH1F("fHistPtBachOmegaPlus", "p_{t} of gen. Bach.;pt_{Bach};Number of Bach.", 200, 0, 10);
     fListHistCascade->Add(fHistPtBachOmegaPlus);
  }
  
  if (!fHistPtMesDghterOmegaPlus) {
     fHistPtMesDghterOmegaPlus = new TH1F("fHistPtMesDghterOmegaPlus", "p_{t} of gen. Meson #Lambda dghter);pt_{MesDght};Number of Mes.", 200, 0, 10);
     fListHistCascade->Add(fHistPtMesDghterOmegaPlus);
  }
    
  if (!fHistPtBarDghterOmegaPlus) {
     fHistPtBarDghterOmegaPlus = new TH1F("fHistPtBarDghterOmegaPlus", "p_{t} of gen. Baryon #Lambda dghter);pt_{BarDght};Number of Bar.", 200, 0, 10);
     fListHistCascade->Add(fHistPtBarDghterOmegaPlus);
  }
    
  
//--------------------------------------------------------------------------------
// Part 2 - Any reconstructed cascades + reconstructed cascades associated with MC
  
		// - Effective mass histos for cascades candidates.
  
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
  
  
  
  		// - Effective mass histos with combined PID
  
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
  
  		// - PID Probability versus MC Pt(bachelor track)
  if(! f2dHistPIDprobaKaonVsMCPtBach ){
	f2dHistPIDprobaKaonVsMCPtBach  = new TH2F( "f2dHistPIDprobaKaonVsMCPtBach" , "Comb. PID proba to be K^{#pm} Vs MC Bach. Pt ; Pt_{MC}(Bach.) (GeV/c); Comb. PID Proba (Bach. = K^{#pm})", 100, 0.0, 5.0, 110, 0.0, 1.10 );
	fListHistCascade->Add(f2dHistPIDprobaKaonVsMCPtBach);
  }
  
  if(! f2dHistPIDprobaPionVsMCPtBach ){
	f2dHistPIDprobaPionVsMCPtBach  = new TH2F( "f2dHistPIDprobaPionVsMCPtBach" , "Comb. PID proba to be #pi^{#pm} Vs MC Bach. Pt ; Pt_{MC}(Bach.) (GeV/c); Comb. PID Proba (Bach. = #pi^{#pm})", 100, 0.0, 5.0, 110, 0.0, 1.10 );
	fListHistCascade->Add(f2dHistPIDprobaPionVsMCPtBach);
  }
  
  
  		// - Effective mass histos with perfect MC PID on the bachelor
  
  if (! fHistMassWithMcPIDXiMinus) {
    fHistMassWithMcPIDXiMinus = new TH1F("fHistMassWithMcPIDXiMinus","#Xi^{-} candidates, with Bach. MC PID;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 400,1.2,2.0);
    fListHistCascade->Add(fHistMassWithMcPIDXiMinus);
  }
  
  if (! fHistMassWithMcPIDXiPlus) {
    fHistMassWithMcPIDXiPlus = new TH1F("fHistMassWithMcPIDXiPlus","#Xi^{+} candidates, with Bach. MC PID;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",400,1.2,2.0);
    fListHistCascade->Add(fHistMassWithMcPIDXiPlus);
  }

  if (! fHistMassWithMcPIDOmegaMinus) {
	fHistMassWithMcPIDOmegaMinus = new TH1F("fHistMassWithMcPIDOmegaMinus","#Omega^{-} candidates, with Bach. MC PID;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassWithMcPIDOmegaMinus);
  }
 
  if (! fHistMassWithMcPIDOmegaPlus) {
	fHistMassWithMcPIDOmegaPlus = new TH1F("fHistMassWithMcPIDOmegaPlus","#Omega^{+} candidates, with Bach. MC PID;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassWithMcPIDOmegaPlus);
  }
  
  
		// - Effective mass histos for cascades candidates ASSOCIATED with MC.
  
  if (! fHistAsMCMassXiMinus) {
	  fHistAsMCMassXiMinus = new TH1F("fHistAsMCMassXiMinus","#Xi^{-} candidates associated to MC;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 400,1.2,2.0);
	  fListHistCascade->Add(fHistAsMCMassXiMinus);
  }
  
  if (! fHistAsMCMassXiPlus) {
	  fHistAsMCMassXiPlus = new TH1F("fHistAsMCMassXiPlus","#Xi^{+} candidates associated to MC;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",400,1.2,2.0);
	  fListHistCascade->Add(fHistAsMCMassXiPlus);
  }

  if (! fHistAsMCMassOmegaMinus) {
	  fHistAsMCMassOmegaMinus = new TH1F("fHistAsMCMassOmegaMinus","#Omega^{-} candidates associated to MC;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500,1.5,2.5);
	  fListHistCascade->Add(fHistAsMCMassOmegaMinus);
  }
 
  if (! fHistAsMCMassOmegaPlus) {
	  fHistAsMCMassOmegaPlus = new TH1F("fHistAsMCMassOmegaPlus","#Omega^{+} candidates associated to MC;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 500,1.5,2.5);
	  fListHistCascade->Add(fHistAsMCMassOmegaPlus);
  }
  
		
		// -  Generated Pt Vs generated Y of the cascade candidates associated with MC 
		//     + having the proper maximum proba of combined PID for the bachelor
  
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
  	f2dHistAsMCGenPtVsGenYXiMinus = new TH2F("f2dHistAsMCGenPtVsGenYXiMinus", "MC P_{t} Vs MC Y of gen. #Xi^{-} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYXiMinus );
  }
  
  if (!f2dHistAsMCGenPtVsGenYXiPlus) {
	  f2dHistAsMCGenPtVsGenYXiPlus = new TH2F("f2dHistAsMCGenPtVsGenYXiPlus", "MC P_{t} Vs MC Y of gen. #Xi^{+} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYXiPlus );
  }
  
  if (!f2dHistAsMCGenPtVsGenYOmegaMinus) {
	  f2dHistAsMCGenPtVsGenYOmegaMinus = new TH2F("f2dHistAsMCGenPtVsGenYOmegaMinus", "MC P_{t} Vs MC Y of gen. #Omega^{-} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYOmegaMinus );
  }
  
  if (!f2dHistAsMCGenPtVsGenYOmegaPlus) {
	  f2dHistAsMCGenPtVsGenYOmegaPlus = new TH2F("f2dHistAsMCGenPtVsGenYOmegaPlus", "MC P_{t} Vs MC Y of gen. #Omega^{+} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYOmegaPlus );
  } 
  
  
  		// - Generated Eta of the the cascade candidates associated with MC
  if (!fHistAsMCGenEtaXiMinus) {
	  fHistAsMCGenEtaXiMinus = new TH1F("fHistAsMCGenEtaXiMinus", "#eta of gen. #Xi^{-} (associated);#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add( fHistAsMCGenEtaXiMinus );
  }
  
  if (!fHistAsMCGenEtaXiPlus) {
	  fHistAsMCGenEtaXiPlus = new TH1F("fHistAsMCGenEtaXiPlus", "#eta of gen. #Xi^{+} (associated);#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add( fHistAsMCGenEtaXiPlus );
  }
  
  if (!fHistAsMCGenEtaOmegaMinus) {
	  fHistAsMCGenEtaOmegaMinus = new TH1F("fHistAsMCGenEtaOmegaMinus", "#eta of gen. #Omega^{-} (associated);#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add( fHistAsMCGenEtaOmegaMinus );
  }
  
  if (!fHistAsMCGenEtaOmegaPlus) {
	  fHistAsMCGenEtaOmegaPlus = new TH1F("fHistAsMCGenEtaOmegaPlus", "#eta of gen. #Omega^{+} (associated);#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add( fHistAsMCGenEtaOmegaPlus );
  }
  
  
  
  		// - Resolution in Pt as function of generated Pt
  
  if(! f2dHistAsMCResPtXiMinus) {
	  f2dHistAsMCResPtXiMinus = new TH2F( "f2dHistAsMCResPtXiMinus", "Resolution in Pt reconstruction for #Xi^{-}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtXiMinus);
  }
  
  if(! f2dHistAsMCResPtXiPlus) {
	  f2dHistAsMCResPtXiPlus = new TH2F( "f2dHistAsMCResPtXiPlus", "Resolution in Pt reconstruction for #Xi^{+}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtXiPlus);
  }
  
  if(! f2dHistAsMCResPtOmegaMinus) {
	  f2dHistAsMCResPtOmegaMinus = new TH2F( "f2dHistAsMCResPtOmegaMinus", "Resolution in Pt reconstruction for #Omega^{-}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtOmegaMinus);
  }
  
  if(! f2dHistAsMCResPtOmegaPlus) {
	  f2dHistAsMCResPtOmegaPlus = new TH2F( "f2dHistAsMCResPtOmegaPlus", "Resolution in Pt reconstruction for #Omega^{+}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtOmegaPlus);
  }
  
  		// - Resolution in R(2D) as function of generated R
  
  if(! f2dHistAsMCResRXiMinus) {
	  f2dHistAsMCResRXiMinus = new TH2F( "f2dHistAsMCResRXiMinus", "Resolution in transv. position for #Xi^{-}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResRXiMinus);
  }
  
  if(! f2dHistAsMCResRXiPlus) {
	  f2dHistAsMCResRXiPlus = new TH2F( "f2dHistAsMCResRXiPlus", "Resolution in transv. position for #Xi^{+}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResRXiPlus);
  }
  
  if(! f2dHistAsMCResROmegaMinus) {
	  f2dHistAsMCResROmegaMinus = new TH2F( "f2dHistAsMCResROmegaMinus", "Resolution in transv. position for #Omega^{-}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResROmegaMinus);
  }
  
  if(! f2dHistAsMCResROmegaPlus) {
	  f2dHistAsMCResROmegaPlus = new TH2F( "f2dHistAsMCResROmegaPlus", "Resolution in transv. position for #Omega^{+}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResROmegaPlus);
  }
  
                // - Resolution in phi as function of generated Pt
    
  if(! f2dHistAsMCResPhiXiMinus) {
          f2dHistAsMCResPhiXiMinus = new TH2F( "f2dHistAsMCResPhiXiMinus", "Resolution in #phi for #Xi^{-}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);
          fListHistCascade->Add(f2dHistAsMCResPhiXiMinus);
  }
  
  if(! f2dHistAsMCResPhiXiPlus) {
          f2dHistAsMCResPhiXiPlus = new TH2F( "f2dHistAsMCResPhiXiPlus", "Resolution in #phi for #Xi^{+}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);
          fListHistCascade->Add(f2dHistAsMCResPhiXiPlus);
  }
  
  if(! f2dHistAsMCResPhiOmegaMinus) {
          f2dHistAsMCResPhiOmegaMinus = new TH2F( "f2dHistAsMCResPhiOmegaMinus", "Resolution in #phi for #Omega^{-}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);  
          fListHistCascade->Add(f2dHistAsMCResPhiOmegaMinus);
  }
  
  if(! f2dHistAsMCResPhiOmegaPlus) {
          f2dHistAsMCResPhiOmegaPlus = new TH2F( "f2dHistAsMCResPhiOmegaPlus", "Resolution in #phi for #Omega^{+}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);
          fListHistCascade->Add(f2dHistAsMCResPhiOmegaPlus);
  }

         //  - Correlation between proton (antiproton) daughter MC pt and Xi/Omega MC pt (to apply Geant/Fluka correction)
    if (!f2dHistAsMCptProtonMCptXiMinus) {
      f2dHistAsMCptProtonMCptXiMinus = new TH2F( "f2dHistAsMCptProtonMCptXiMinus", "Proton MC pt vs Xi- MC pt", 100, 0., 10., 100, 0., 10.); 
      fListHistCascade->Add(f2dHistAsMCptProtonMCptXiMinus);
    }
    if (!f2dHistAsMCptAntiprotonMCptXiPlus) {
      f2dHistAsMCptAntiprotonMCptXiPlus = new TH2F( "f2dHistAsMCptAntiprotonMCptXiPlus", "Antiproton MC pt vs Xi+ MC pt", 100, 0., 10., 100, 0., 10.);
      fListHistCascade->Add(f2dHistAsMCptAntiprotonMCptXiPlus);
    }
    if (!f2dHistAsMCptProtonMCptOmegaMinus) {
      f2dHistAsMCptProtonMCptOmegaMinus = new TH2F( "f2dHistAsMCptProtonMCptOmegaMinus", "Proton MC pt vs Omega- MC pt", 100, 0., 10., 100, 0., 10.);
      fListHistCascade->Add(f2dHistAsMCptProtonMCptOmegaMinus);
    }
    if (!f2dHistAsMCptAntiprotonMCptOmegaPlus) {
      f2dHistAsMCptAntiprotonMCptOmegaPlus = new TH2F( "f2dHistAsMCptAntiprotonMCptOmegaPlus", "Antiproton MC pt vs Omega+ MC pt", 100, 0., 10., 100, 0., 10.);
      fListHistCascade->Add(f2dHistAsMCptAntiprotonMCptOmegaPlus);
    }

    if (! fHistV0toXiCosineOfPointingAngle) {
        fHistV0toXiCosineOfPointingAngle = new TH1F("fHistV0toXiCosineOfPointingAngle", "Cos. of V0 Ptng Angl / Xi vtx ;Cos(V0 Point. Angl / Xi vtx); Counts", 200, 0.95, 1.0001);
        fListHistCascade->Add(fHistV0toXiCosineOfPointingAngle);
    }

    if (! fHistV0CosineOfPointingAnglevsPtXi) {
        fHistV0CosineOfPointingAnglevsPtXi = new TH2F("fHistV0CosineOfPointingAnglevsPtXi", "Cos. of V0 Ptng Angl vs cascade Pt ;Cascade p_{T} (GeV/c); Cos(V0PA) wrt Xi vtx",100, 0., 10., 251, 0.95, 1.0002);
        fListHistCascade->Add(fHistV0CosineOfPointingAnglevsPtXi);
    }
    if (! fHistV0CosineOfPointingAnglevsPtOmega) {
        fHistV0CosineOfPointingAnglevsPtOmega = new TH2F("fHistV0CosineOfPointingAnglevsPtOmega", "Cos. of V0 Ptng Angl vs cascade Pt ;Cascade p_{T} (GeV/c); Cos(V0PA) wrt Omega vtx",100, 0., 10., 251, 0.95, 1.0002);
        fListHistCascade->Add(fHistV0CosineOfPointingAnglevsPtOmega);
    }
                
    // - PID container
if(! fCFContCascadePIDAsXiMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 800;
  lNbBinsPerVar[2] = 22;
  lNbBinsPerVar[3] = 11;

  fCFContCascadePIDAsXiMinus = new AliCFContainer("fCFContCascadePIDAsXiMinus","Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  //setting the bin limits 
  fCFContCascadePIDAsXiMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDAsXiMinus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
  fCFContCascadePIDAsXiMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  Double_t *lBinLim3  = new Double_t[ lNbBinsPerVar[3]+1 ];
  for(Int_t i=3; i< lNbBinsPerVar[3]+1;i++)   lBinLim3[i]  = (Double_t)(i-1)*10.;
  lBinLim3[0] = 0.0;
  lBinLim3[1] = 5.0;
  lBinLim3[2] = 10.0;
  fCFContCascadePIDAsXiMinus->SetBinLimits(3,  lBinLim3 );       // Centrality
  
  // Setting the step title : one per PID case
  fCFContCascadePIDAsXiMinus->SetStepTitle(0, "No PID");
  fCFContCascadePIDAsXiMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDAsXiMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDAsXiMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDAsXiMinus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDAsXiMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDAsXiMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDAsXiMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDAsXiMinus->SetVarTitle(1, "M( #Lambda , #pi^{-} ) (GeV/c^{2})");
  fCFContCascadePIDAsXiMinus->SetVarTitle(2, "Y_{#Xi}");
  fCFContCascadePIDAsXiMinus->SetVarTitle(3, "Centrality");
}

if(! fCFContCascadePIDAsXiPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 800;
  lNbBinsPerVar[2] = 22;
  lNbBinsPerVar[3] = 11;
  
  fCFContCascadePIDAsXiPlus = new AliCFContainer("fCFContCascadePIDAsXiPlus","Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDAsXiPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDAsXiPlus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
  fCFContCascadePIDAsXiPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  Double_t *lBinLim3  = new Double_t[ lNbBinsPerVar[3]+1 ];
  for(Int_t i=3; i< lNbBinsPerVar[3]+1;i++)   lBinLim3[i]  = (Double_t)(i-1)*10.;
  lBinLim3[0] = 0.0;
  lBinLim3[1] = 5.0;
  lBinLim3[2] = 10.0;
  fCFContCascadePIDAsXiPlus->SetBinLimits(3,lBinLim3);     // Centrality
  
  // Setting the step title : one per PID case
  fCFContCascadePIDAsXiPlus->SetStepTitle(0, "No PID");
  fCFContCascadePIDAsXiPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDAsXiPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDAsXiPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDAsXiPlus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDAsXiPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDAsXiPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDAsXiPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDAsXiPlus->SetVarTitle(1, "M( #Lambda , #pi^{+} ) (GeV/c^{2})");
  fCFContCascadePIDAsXiPlus->SetVarTitle(2, "Y_{#Xi}");
  fCFContCascadePIDAsXiPlus->SetVarTitle(3, "Centrality");
}


if(! fCFContCascadePIDAsOmegaMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 1000;
  lNbBinsPerVar[2] = 22;
  lNbBinsPerVar[3] = 11;
  
  fCFContCascadePIDAsOmegaMinus = new AliCFContainer("fCFContCascadePIDAsOmegaMinus","Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits 
  fCFContCascadePIDAsOmegaMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDAsOmegaMinus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
  fCFContCascadePIDAsOmegaMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  Double_t *lBinLim3  = new Double_t[ lNbBinsPerVar[3]+1 ];
  for(Int_t i=3; i< lNbBinsPerVar[3]+1;i++)   lBinLim3[i]  = (Double_t)(i-1)*10.;
  lBinLim3[0] = 0.0;
  lBinLim3[1] = 5.0;
  lBinLim3[2] = 10.0;
  fCFContCascadePIDAsOmegaMinus->SetBinLimits(3,lBinLim3);     // Centrality
  
  // Setting the step title : one per PID case
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(0, "No PID");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDAsOmegaMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDAsOmegaMinus->SetVarTitle(1, "M( #Lambda , K^{-} ) (GeV/c^{2})");
  fCFContCascadePIDAsOmegaMinus->SetVarTitle(2, "Y_{#Omega}");
  fCFContCascadePIDAsOmegaMinus->SetVarTitle(3, "Centrality");
}

if(! fCFContCascadePIDAsOmegaPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4]= {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 1000;
  lNbBinsPerVar[2] = 22;
  lNbBinsPerVar[3] = 11;
  
  fCFContCascadePIDAsOmegaPlus = new AliCFContainer("fCFContCascadePIDAsOmegaPlus","Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits 
  fCFContCascadePIDAsOmegaPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDAsOmegaPlus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
  fCFContCascadePIDAsOmegaPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  Double_t *lBinLim3  = new Double_t[ lNbBinsPerVar[3]+1 ];
  for(Int_t i=3; i< lNbBinsPerVar[3]+1;i++)   lBinLim3[i]  = (Double_t)(i-1)*10.;
  lBinLim3[0] = 0.0;
  lBinLim3[1] = 5.0;
  lBinLim3[2] = 10.0;
  fCFContCascadePIDAsOmegaPlus->SetBinLimits(3,lBinLim3);     // Centrality
  
  // Setting the step title : one per PID case
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(0, "No PID");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDAsOmegaPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDAsOmegaPlus->SetVarTitle(1, "M( #Lambda , K^{+} ) (GeV/c^{2})");
  fCFContCascadePIDAsOmegaPlus->SetVarTitle(2, "Y_{#Omega}");
  fCFContCascadePIDAsOmegaPlus->SetVarTitle(3, "Centrality");
  
//  fListHistCascade->Add(fCFContCascadePIDAsOmegaPlus);
  
}

// Part 3 : Towards the optimisation of topological selections -------
if(! fCFContAsCascadeCuts){
   
	// Container meant to store all the relevant distributions corresponding to the cut variables.
        //          - NB overflow/underflow of variables on which we want to cut later should be 0!!!

  const	Int_t  lNbSteps      =  4 ;
  const Int_t  lNbVariables  =  22 ;
  
  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[lNbVariables] = {0};
  lNbBinsPerVar[0]  = 100;
  lNbBinsPerVar[1]  = 126;
  lNbBinsPerVar[2]  = 24;
  lNbBinsPerVar[3]  = 220;
  lNbBinsPerVar[4]  = 30;
  lNbBinsPerVar[5]  = 50;
  
  lNbBinsPerVar[6]  = 101;
  lNbBinsPerVar[7]  = 102;
  lNbBinsPerVar[8]  = 101;
  lNbBinsPerVar[9]  = 26;
  lNbBinsPerVar[10] = 26;
  
  lNbBinsPerVar[11] = 150; // 1-MeV/c2 bins
  lNbBinsPerVar[12] = 120;  
  
  lNbBinsPerVar[13] = 100;
  
  lNbBinsPerVar[14] = 110; // 0.02 unit of y per bin 
  lNbBinsPerVar[15] = 110; 
  
  lNbBinsPerVar[16] = 2;
 
  lNbBinsPerVar[17] = 11;
  lNbBinsPerVar[18] = 100;

  lNbBinsPerVar[19] = 112; // Proper time of cascade
  lNbBinsPerVar[20] = 112; // Proper time of V0 
  lNbBinsPerVar[21] = 112; // Distance V0-Xi in the transverse plane 
  
  fCFContAsCascadeCuts = new AliCFContainer("fCFContAsCascadeCuts","Cut Container for Asso. Cascades", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  //0
  fCFContAsCascadeCuts->SetBinLimits(0, 0., 2.);                 // DcaXiDaughters : 0.0 to 2.0
  //1
   Double_t *lBinLim1  = new Double_t[ lNbBinsPerVar[1]+1 ];
   lBinLim1[0] = 0.0;
   lBinLim1[1] = 0.03;
   for(Int_t i=2; i< lNbBinsPerVar[1];i++)   lBinLim1[i]  = (Double_t)0.03   + (5.  - 0.03 )/(lNbBinsPerVar[1]-2)  * (Double_t)(i-1) ;
   lBinLim1[ lNbBinsPerVar[1]  ] = 100.0;
   fCFContAsCascadeCuts -> SetBinLimits(1,  lBinLim1 );
  delete [] lBinLim1;                                            // DcaBachToPrimVertexXi : 0.0 to 0.5
  //2
  fCFContAsCascadeCuts->SetBinLimits(2, .9988, 1.);              // XiCosineOfPointingAngle : 0.99 to 1.0        
  //3
  fCFContAsCascadeCuts -> SetBinLimits(3, 0., 110. );            // XiRadius : 0.0 to 110.0

  //4
  fCFContAsCascadeCuts->SetBinLimits(4, 1.1, 1.13);              // InvMassLambdaAsCascDghter
  //5
  fCFContAsCascadeCuts->SetBinLimits(5, 0., 2.);                 // DcaV0DaughtersXi : 0.0 to 2.0        
  //6
  fCFContAsCascadeCuts->SetBinLimits(6, 0.98, 1.0002);           // V0CosineOfPointingAngleXi : 0.98 to 1.0     
  //7
  Double_t *lBinLim7  = new Double_t[ lNbBinsPerVar[7]+1 ];
  for(Int_t i=0; i< lNbBinsPerVar[7]-1;i++)   lBinLim7[i]  = (Double_t)0.0   + (100.  - 0.0 )/(lNbBinsPerVar[7]-2)  * (Double_t)i ;
  lBinLim7[ lNbBinsPerVar[7]-1] = 200.0;
  lBinLim7[ lNbBinsPerVar[7]] = 1000.0;
  fCFContAsCascadeCuts -> SetBinLimits(7,  lBinLim7 );
  delete [] lBinLim7;                                             // V0RadiusXi : 0.0 to 1000.0      
  //8
  Double_t *lBinLim8  = new Double_t[ lNbBinsPerVar[8]+1 ];
  for(Int_t i=0; i< lNbBinsPerVar[8];i++)   lBinLim8[i]  = (Double_t)0.0   + (0.4  - 0.0 )/(lNbBinsPerVar[8]-1)  * (Double_t)i ;
  lBinLim8[ lNbBinsPerVar[8]  ] = 100.0;
  fCFContAsCascadeCuts -> SetBinLimits(8,  lBinLim8 );
  delete [] lBinLim8;                                             // DcaV0ToPrimVertexXi : 0. to 0.4     
  //9
  Double_t *lBinLim9  = new Double_t[ lNbBinsPerVar[9]+1 ];
  for(Int_t i=0; i< lNbBinsPerVar[9];i++)   lBinLim9[i]  = (Double_t)0.0   + (0.25  - 0.0 )/(lNbBinsPerVar[9]-1)  * (Double_t)i ;
  lBinLim9[ lNbBinsPerVar[9]  ] = 100.0;
  fCFContAsCascadeCuts -> SetBinLimits(9,  lBinLim9 );
  delete [] lBinLim9;                                             // DcaPosToPrimVertexXi : 0.0 to 0.25   
  //10
  Double_t *lBinLim10  = new Double_t[ lNbBinsPerVar[10]+1 ];
  for(Int_t i=0; i< lNbBinsPerVar[10];i++)   lBinLim10[i]  = (Double_t)0.0   + (0.25  - 0.0 )/(lNbBinsPerVar[10]-1)  * (Double_t)i ;
  lBinLim10[ lNbBinsPerVar[10]  ] = 100.0;
  fCFContAsCascadeCuts -> SetBinLimits(10,  lBinLim10 );
  delete [] lBinLim10;                                            // DcaPosToPrimVertexXi : 0.0 to 0.25 
  //11
  fCFContAsCascadeCuts->SetBinLimits(11, 1.25, 1.40);	          // InvMassXi
  fCFContAsCascadeCuts->SetBinLimits(12, 1.62, 1.74);	          // InvMassOmega
  fCFContAsCascadeCuts->SetBinLimits(13, 0.0, 10.0);	          // XiTransvMom 
  fCFContAsCascadeCuts->SetBinLimits(14, -1.1, 1.1);              // Y(Xi)
  fCFContAsCascadeCuts->SetBinLimits(15, -1.1, 1.1);              // Y(Omega)
  fCFContAsCascadeCuts->SetBinLimits(16, 0.0, 2.0);               // Natural or injected  
  Double_t *lBinLim17  = new Double_t[ lNbBinsPerVar[17]+1 ];
  for(Int_t i=3; i< lNbBinsPerVar[17]+1;i++)   lBinLim17[i]  = (Double_t)(i-1)*10.;
  lBinLim17[0] = 0.0;
  lBinLim17[1] = 5.0;
  lBinLim17[2] = 10.0;
  fCFContAsCascadeCuts -> SetBinLimits(17,  lBinLim17 );       // Centrality
  delete [] lBinLim17;
  fCFContAsCascadeCuts->SetBinLimits(18, 0.0, 6000.0);         // ESD track multiplicity 
  Double_t *lBinLim19  = new Double_t[ lNbBinsPerVar[19]+1 ];
  for(Int_t i=0; i< lNbBinsPerVar[19];i++)   lBinLim19[i]  = (Double_t)-1.   + (110.  + 1.0 )/(lNbBinsPerVar[19]-1)  * (Double_t)i ;
  lBinLim19[ lNbBinsPerVar[19]  ] = 2000.0;
  fCFContAsCascadeCuts->SetBinLimits(19, lBinLim19);           // Proper time cascade 
  fCFContAsCascadeCuts->SetBinLimits(20, lBinLim19);           // Proper time V0
  fCFContAsCascadeCuts->SetBinLimits(21, lBinLim19);           // Distance V0-Xi in the transverse plane
  delete [] lBinLim19;
  // Setting the number of steps : one for each cascade species (Xi-, Xi+ and Omega-, Omega+)
  fCFContAsCascadeCuts->SetStepTitle(0, "#Xi^{-} candidates associated to MC");
  fCFContAsCascadeCuts->SetStepTitle(1, "#bar{#Xi}^{+} candidates associated to MC");
  fCFContAsCascadeCuts->SetStepTitle(2, "#Omega^{-} candidates associated to MC");
  fCFContAsCascadeCuts->SetStepTitle(3, "#bar{#Omega}^{+} candidates associated to MC");
  
  // Setting the variable title, per axis
  fCFContAsCascadeCuts->SetVarTitle(0,  "DCA(XiDaughters) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(1,  "DCA(Bach/PrimVertex) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(2,  "cos(Xi pointing angle)");
  fCFContAsCascadeCuts->SetVarTitle(3,  "R_{2d}(Xi decay) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(4,  "M_{#Lambda}(As Casc Dghter) (GeV/c^{2})");
  fCFContAsCascadeCuts->SetVarTitle(5,  "DCA(V0 Daughters Xi) (cm)");
  
  fCFContAsCascadeCuts->SetVarTitle(6,  "cos(V0 pointing Angle) in Casc");
  fCFContAsCascadeCuts->SetVarTitle(7,  "R_{2d}(V0 decay) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(8,  "DCA(V0/PrimVertex) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(9,  "DCA(Pos/PrimVertex) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(10, "DCA(Neg/PrimVertex) (cm)");
  
  fCFContAsCascadeCuts->SetVarTitle(11, "Inv. Mass(Xi) (GeV/c^{2})");
  fCFContAsCascadeCuts->SetVarTitle(12, "Inv. Mass(Omega) (GeV/c^{2})");
  
  fCFContAsCascadeCuts->SetVarTitle(13, "Pt_{MC}(Casc.) (GeV/c)");
  
  fCFContAsCascadeCuts->SetVarTitle(14, "Y_{MC}(Xi)");
  fCFContAsCascadeCuts->SetVarTitle(15, "Y_{MC}(Omega)");
  
  fCFContAsCascadeCuts->SetVarTitle(16, "Injected or natural particles");
  
  fCFContAsCascadeCuts->SetVarTitle(17, "Centrality");
  fCFContAsCascadeCuts->SetVarTitle(18, "ESD track multiplicity");
  
  fCFContAsCascadeCuts->SetVarTitle(19, "Proper time cascade");
  fCFContAsCascadeCuts->SetVarTitle(20, "Proper time V0"); 
  fCFContAsCascadeCuts->SetVarTitle(21, "Distance V0-Xi in the transverse plane");
}

  fV0Ampl = new TH1F("fV0Ampl","",500,0.,30000);
  fListHistCascade->Add(fV0Ampl);

  fHistEtaBachXiM = new TH1F("fHistEtaBachXiM","",40,-2.,2.);
  fListHistCascade->Add(fHistEtaBachXiM);
  fHistEtaPosXiM = new TH1F("fHistEtaPosXiM","",40,-2.,2.);
  fListHistCascade->Add(fHistEtaPosXiM); 
  fHistEtaNegXiM = new TH1F("fHistEtaNegXiM","",40,-2.,2.);
  fListHistCascade->Add(fHistEtaNegXiM); 
 


PostData(1, fListHistCascade); 
PostData(2, fCFContCascadePIDAsXiMinus);
PostData(3, fCFContCascadePIDAsXiPlus);
PostData(4, fCFContCascadePIDAsOmegaMinus);
PostData(5, fCFContCascadePIDAsOmegaPlus);
PostData(6, fCFContAsCascadeCuts);

}// end CreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskCheckPerformanceCascadePbPb::UserExec(Option_t *) {
	
  // Main loop
  // Called for each event
	
	AliESDEvent *lESDevent = 0x0;
	AliAODEvent *lAODevent = 0x0;
	AliMCEvent  *lMCevent  = 0x0; 
	AliStack    *lMCstack  = 0x0; 

        TClonesArray *arrayMC = 0;

        if (!fPIDResponse) {
          AliError("Cannot get pid response");
          return;
        }

	
  // Connect to the InputEvent	
  // After these lines, we should have an ESD/AOD event 
		
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

          if (fkRerunV0CascVertexers) { // relaunch V0 and Cascade vertexer
              cout<<" ---- Re-vertexing ---- "<<endl;
              lESDevent->ResetCascades();
              lESDevent->ResetV0s();
              AliV0vertexer lV0vtxer;
              AliCascadeVertexer lCascVtxer;
              lV0vtxer.SetCuts(fV0VertexerSels);
              lCascVtxer.SetCuts(fCascadeVertexerSels);
              cout<<"  - V0 #chi^{2} _________________ <"<<fV0VertexerSels[0]<<endl;
              cout<<"  - DCA(prim. Vtx/ 1^{st} daughter) ___ >"<<fV0VertexerSels[1]<<endl;
              cout<<"  - DCA(prim. Vtx/ 2^{nd} daughter) __  >"<<fV0VertexerSels[2]<<endl;
              cout<<"  - DCA between V0 daughters ___ <"<<fV0VertexerSels[3]<<endl;
              cout<<"  - cos(V0 pointing angle) _______ >"<<fV0VertexerSels[4]<<endl;
              cout<<"  - R_{transv}(V0 decay) ________ >"<<fV0VertexerSels[5]<<endl;
              cout<<"  - R_{transv}(V0 decay) ________ <"<<fV0VertexerSels[6]<<endl;
              cout<<"  - Casc. #chi^{2} ______________  <"<<fCascadeVertexerSels[0]<<endl;
              cout<<"  - DCA(prim. Vtx/ V0) _________ >"<<fCascadeVertexerSels[1]<<endl;
              cout<<"  - | M_{#Lambda}(reco) - M_{#Lambda}(pdg) | _______ <"<<fCascadeVertexerSels[2]<<endl;
              cout<<"  - DCA(prim. Vtx/ Bach) _______ >"<<fCascadeVertexerSels[3]<<endl;
              cout<<"  - DCA between Bach/ #Lambda ______ <"<<fCascadeVertexerSels[4]<<endl;
              cout<<"  - cos(Casc. pointing angle) ____ >"<<fCascadeVertexerSels[5]<<endl;
              cout<<"  - R_{transv}(Casc. decay) ______ >"<<fCascadeVertexerSels[6]<<endl;
              cout<<"  - R_{transv}(Casc. decay) ______ <"<<fCascadeVertexerSels[7]<<endl;
              lV0vtxer.Tracks2V0vertices(lESDevent);
              lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
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
	

  //------------------------------------------------
  // 1 - Preparing the general info about of the event = prim. Vtx + magnetic field (ESD)
 
        Double_t lBestPrimaryVtxPos[3]   = {-100.0, -100.0, -100.0};
        Double_t lMagneticField    = -10.; 
        Int_t ncascades = -1;
        // Primary tracks from ESD/AOD
        Float_t lPrimaryTrackMultiplicity= -1.;

//        Int_t    lSPDTrackletsMultiplicity = -1;

        AliCentrality* centrality = 0x0;
        AliESDVZERO* esdV0 = 0x0;
        AliAODVZERO* aodV0 = 0x0;
        Float_t multV0A = 0.;
        Float_t multV0C = 0.;

        Int_t   nNumberOfMCPrimaries     = -1;

        Int_t nTrackMultiplicity = (InputEvent())->GetNumberOfTracks();

        if (fAnalysisType == "ESD" ) {

        // Prim vertex
          const AliESDVertex *lPrimaryTrackingVtx = lESDevent->GetPrimaryVertexTracks();  // get the vtx stored in ESD found with tracks
          const AliESDVertex *lPrimarySPDVtx      = lESDevent->GetPrimaryVertexSPD();     // get the vtx stored in ESD found with SPD tracklets

          const AliESDVertex *lPrimaryBestVtx     = lESDevent->GetPrimaryVertex();
          // get the best primary vertex available for the event
          // As done in AliCascadeVertexer, we keep the one which is the best one available.
          // between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
          lPrimaryBestVtx->GetXYZ( lBestPrimaryVtxPos );

          // FIXME : quality selection regarding pile-up rejection, not in the Cascade task 
          if (fkRejectEventPileUp) {
            if (lESDevent->IsPileupFromSPD() ) {// minContributors=3, minZdist=0.8, nSigmaZdist=3., nSigmaDiamXY=2., nSigmaDiamZ=5.  -> see http://alisoft.cern.ch/viewvc/trunk/STEER/AliESDEvent.h?root=AliRoot&r1=41914&r2=42199&pathrev=42199
              AliWarning("Pb / Event tagged as pile-up by SPD... return !"); 
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDAsXiMinus);
              PostData(3, fCFContCascadePIDAsXiPlus);
              PostData(4, fCFContCascadePIDAsOmegaMinus);
              PostData(5, fCFContCascadePIDAsOmegaPlus);
              PostData(6, fCFContAsCascadeCuts);
 
              return; 
            }
          }
          // FIXME : remove TPC-only primary vertex : retain only events with tracking + SPD vertex
          if (fkQualityCutNoTPConlyPrimVtx) {
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

          lMagneticField = lESDevent->GetMagneticField();

          ncascades = lESDevent->GetNumberOfCascades();

          lPrimaryTrackMultiplicity = fESDtrackCuts->CountAcceptedTracks(lESDevent);

//          const AliMultiplicity *lAliMult = lESDevent->GetMultiplicity();
//          lSPDTrackletsMultiplicity       = lAliMult->GetNumberOfTracklets();

          //        nNumberOfMCPrimaries       = lMCstack->GetNprimary();
          nNumberOfMCPrimaries = lMCstack->GetNtrack(); // MN: this stack->GetNtrack(); has to be used to check primaries because in HIJING decay products of D and B mesons are also primaries and produced in HIJING during transport, however this is not the number of primaries!
          // Centrality determination
          esdV0 = lESDevent->GetVZEROData();
          multV0A = esdV0->GetMTotV0A();
          multV0C = esdV0->GetMTotV0C();

          centrality = lESDevent->GetCentrality();

        } else if (fAnalysisType == "AOD") {

          const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();
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

          // get the best primary vertex available for the event GetVertex(0)
          // This one will be used for next calculations (DCA essentially)
          lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );

          lMagneticField = lAODevent->GetMagneticField();

          ncascades = lAODevent->GetNumberOfCascades();

          lPrimaryTrackMultiplicity = 0.;
          for (Int_t itrack = 0; itrack<nTrackMultiplicity; itrack++) {
            AliAODTrack* track = dynamic_cast<AliAODTrack*>(lAODevent->GetTrack(itrack));
            if(!track) AliFatal("Not a standard AOD");
            if (track->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) lPrimaryTrackMultiplicity++;
          }

//          lSPDTrackletsMultiplicity = lAODevent->GetTracklets()->GetNumberOfTracklets();

          nNumberOfMCPrimaries = arrayMC->GetEntries();


          centrality = lAODevent->GetCentrality();
          aodV0 = lAODevent->GetVZEROData();
          multV0A=aodV0->GetMTotV0A();
          multV0C=aodV0->GetMTotV0C();

        } 

        if (fkQualityCutZprimVtxPos) { // Should be already in the centrality selection
          if (TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRange ) {
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
 

      //  Printf("Centrality percentile V0M for this event %f)\n",  centrality->GetCentralityPercentile("V0M"));


        Float_t lcentrality = 0.;
        
        if (fkUseCleaning) lcentrality = centrality->GetCentralityPercentile(fCentrEstimator.Data());
        else {
          lcentrality = centrality->GetCentralityPercentileUnchecked(fCentrEstimator.Data());
          if (centrality->GetQuality()>1) {
            PostData(1, fListHistCascade);
            PostData(2, fCFContCascadePIDAsXiMinus);
            PostData(3, fCFContCascadePIDAsXiPlus);
            PostData(4, fCFContCascadePIDAsOmegaMinus);
            PostData(5, fCFContCascadePIDAsOmegaPlus);
            PostData(6, fCFContAsCascadeCuts);

            return;
          } 
        }
        if (lcentrality<fCentrLowLim||lcentrality>=fCentrUpLim) {
          PostData(1, fListHistCascade);
          return;
        }

        if (nNumberOfMCPrimaries < 1) {
          PostData(1, fListHistCascade);
          PostData(2, fCFContCascadePIDAsXiMinus);
          PostData(3, fCFContCascadePIDAsXiPlus);
          PostData(4, fCFContCascadePIDAsOmegaMinus);
          PostData(5, fCFContCascadePIDAsOmegaPlus);
          PostData(6, fCFContAsCascadeCuts);
  
          return; // should be useless because we require vertex and centrality selection 
        }

        fV0Ampl->Fill(multV0A+multV0C);

        // Best vertex distribution
        fHistBestVtxX->Fill(lBestPrimaryVtxPos[0]);
        fHistBestVtxY->Fill(lBestPrimaryVtxPos[1]);
        fHistBestVtxZ->Fill(lBestPrimaryVtxPos[2]);
 
        fHistEvtsInCentralityBinsvsNtracks->Fill(lcentrality,lPrimaryTrackMultiplicity);

        fHistMCTrackMultiplicity->Fill( nNumberOfMCPrimaries );  // MN: neutral particles included and also not physical ones

        //cout << "Name of the accessed file :" <<  CurrentFileName() << endl;

        Int_t   nMCPrimariesInAcceptance   =  0;
//_____________________________________________________________________________	
// Part 1 - Loop over the MC primaries	
        
    for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < nNumberOfMCPrimaries; iCurrentLabelStack++) {

      if (fAnalysisType == "ESD") {

        TParticle* lCurrentParticle = 0x0; 
                   lCurrentParticle = lMCstack->Particle( iCurrentLabelStack );
        if (!lCurrentParticle) {
          Printf("MC Primary loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
          continue;
        }
  
        if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) continue;
        
       
        if ( TMath::Abs( lCurrentParticle->Eta() ) > 0.8 ) continue;    

      } else if (fAnalysisType == "AOD") {

        AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack);
        if(TMath::Abs(partMC->Eta())>0.8) continue;  

      }

      nMCPrimariesInAcceptance++;
 
   }
    
   f2dHistRecoMultVsMCMult->Fill( lPrimaryTrackMultiplicity, nMCPrimariesInAcceptance );  // MN: neutral are included


   // For proton
   
/*   for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < nNumberOfMCPrimaries; iCurrentLabelStack++) {
          
    	TParticle* lCurrentParticle = lMCstack->Particle( iCurrentLabelStack );
	if(!lCurrentParticle){
		Printf("Proton loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
		continue;
		
	}
	
	if( lCurrentParticle->GetPdgCode() == 2212 )
		fHistEtaGenProton->Fill( lCurrentParticle->Eta() );

	if( lCurrentParticle->GetPdgCode() == -2212 )
		fHistEtaGenAntiProton->Fill( lCurrentParticle->Eta() );
      }// end loop over primary proton
*/   

      
       
//_____________________________________________________________________________	
// Part 2 - Loop over the different types of GENERATED cascades (Xi-+, Omega-+)	

	// - Initialisation of useful local variables
		
	Int_t lPdgCodeCasc            = 0;
	Int_t lPdgCodeBach            = 0;
	Int_t lPdgCodeLambda          = 0;
	Int_t lPdgCodeDghtMesV0       = 0;
	Int_t lPdgCodeDghtBarV0       = 0;
	
	
	TH1F *lHistEtaGenCasc         = 0;	
	TH3D *l3dHistGenPtVsGenYvsCentNat = 0;
        TH3D *l3dHistGenPtVsGenYvsNtracksNat = 0;
        TH3D *l3dHistGenPtVsGenYvsCentInj = 0;
        TH3D *l3dHistGenPtVsGenYvsNtracksInj = 0;
        TH3D *l3dHistGenPtVsGenctauvsCentNat = 0;
        TH3D *l3dHistGenPtVsGenctauvsCentInj = 0; 

	TH1F *lHistThetaGenCascNat    = 0;
        TH1F *lHistThetaGenCascInj    = 0;
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

        Int_t endOfHijingEvent = -1;

for (Int_t iCascType = 1; iCascType < 5; iCascType++) { 
  ncascperev = 0;
  ncascperevtot = 0;

  endOfHijingEvent = -1;
       
  switch (iCascType) {
    case 1: // Xi-
         lPdgCodeCasc       =   3312;  //Xi-
         lPdgCodeBach       =   -211;  //Pi-
         lPdgCodeLambda     =   3122;  //Lambda0
         lPdgCodeDghtMesV0  =   -211;  //Pi-
         lPdgCodeDghtBarV0  =   2212;  //Proton 
	 	
	 	// any Xi-
	 lHistEtaGenCasc        = fHistEtaGenCascXiMinus;
         l3dHistGenPtVsGenYvsCentNat = f3dHistGenPtVsGenYvsCentXiMinusNat;
         l3dHistGenPtVsGenYvsNtracksNat = f3dHistGenPtVsGenYvsNtracksXiMinusNat;
         l3dHistGenPtVsGenYvsCentInj = f3dHistGenPtVsGenYvsCentXiMinusInj;
         l3dHistGenPtVsGenYvsNtracksInj = f3dHistGenPtVsGenYvsNtracksXiMinusInj;
         l3dHistGenPtVsGenctauvsCentNat = f3dHistGenPtVsGenctauvsCentXiMinusNat;
         l3dHistGenPtVsGenctauvsCentInj = f3dHistGenPtVsGenctauvsCentXiMinusInj;

	 	// cascades generated within acceptance (cut in pt + theta)
	 lHistThetaGenCascNat      = fHistThetaGenCascXiMinusNat;
         lHistThetaGenCascInj      = fHistThetaGenCascXiMinusInj;
	 l2dHistGenPtVsGenYFdbl = f2dHistGenPtVsGenYFdblXiMinus;
	 lHistThetaLambda       = fHistThetaLambdaXiMinus;
	 lHistThetaBach         = fHistThetaBachXiMinus;
	 lHistThetaBarDghter    = fHistThetaBarDghterXiMinus;
	 lHistThetaMesDghter    = fHistThetaMesDghterXiMinus;
	 lHistPtBach	        = fHistPtBachXiMinus;
	 lHistPtBarDghter       = fHistPtBarDghterXiMinus;
	 lHistPtMesDghter       = fHistPtMesDghterXiMinus;
        break; 
           
    case 2: // Xi+
         lPdgCodeCasc        =  -3312;  //Xi+
         lPdgCodeBach        =    211;  //Pi+
         lPdgCodeLambda      =  -3122;  //AntiLambda0
         lPdgCodeDghtMesV0   =    211;  //Pi+
         lPdgCodeDghtBarV0   =  -2212;  //AntiProton  
	 
	 	// any Xi+
	 lHistEtaGenCasc     	= fHistEtaGenCascXiPlus;
         l3dHistGenPtVsGenYvsCentNat = f3dHistGenPtVsGenYvsCentXiPlusNat;
         l3dHistGenPtVsGenYvsNtracksNat = f3dHistGenPtVsGenYvsNtracksXiPlusNat;
         l3dHistGenPtVsGenYvsCentInj = f3dHistGenPtVsGenYvsCentXiPlusInj;
         l3dHistGenPtVsGenYvsNtracksInj = f3dHistGenPtVsGenYvsNtracksXiPlusInj;
         l3dHistGenPtVsGenctauvsCentNat = f3dHistGenPtVsGenctauvsCentXiPlusNat;
         l3dHistGenPtVsGenctauvsCentInj = f3dHistGenPtVsGenctauvsCentXiPlusInj;
	
	 	// cascades generated within acceptance (cut in pt + theta)	 
	 lHistThetaGenCascNat      = fHistThetaGenCascXiPlusNat;
         lHistThetaGenCascInj      = fHistThetaGenCascXiPlusInj;
	 l2dHistGenPtVsGenYFdbl = f2dHistGenPtVsGenYFdblXiPlus;
	 lHistThetaLambda       = fHistThetaLambdaXiPlus;
	 lHistThetaBach         = fHistThetaBachXiPlus;
	 lHistThetaBarDghter    = fHistThetaBarDghterXiPlus;
	 lHistThetaMesDghter    = fHistThetaMesDghterXiPlus;
	 lHistPtBach	        = fHistPtBachXiPlus;
	 lHistPtBarDghter       = fHistPtBarDghterXiPlus;
	 lHistPtMesDghter       = fHistPtMesDghterXiPlus;  
    	break;
   
    case 3: // Omega-
    	 lPdgCodeCasc       =   3334;  //Omega-
         lPdgCodeBach       =   -321;  //K-
         lPdgCodeLambda     =   3122;  //Lambda0
         lPdgCodeDghtMesV0  =   -211;  //Pi-
         lPdgCodeDghtBarV0  =   2212;  //Proton 
	 
	 	// any Omega-
	 lHistEtaGenCasc        = fHistEtaGenCascOmegaMinus;	 	
         l3dHistGenPtVsGenYvsCentNat = f3dHistGenPtVsGenYvsCentOmegaMinusNat;
         l3dHistGenPtVsGenYvsNtracksNat = f3dHistGenPtVsGenYvsNtracksOmegaMinusNat;
         l3dHistGenPtVsGenYvsCentInj = f3dHistGenPtVsGenYvsCentOmegaMinusInj;
         l3dHistGenPtVsGenYvsNtracksInj = f3dHistGenPtVsGenYvsNtracksOmegaMinusInj; 
         l3dHistGenPtVsGenctauvsCentNat = f3dHistGenPtVsGenctauvsCentOmegaMinusNat;
         l3dHistGenPtVsGenctauvsCentInj = f3dHistGenPtVsGenctauvsCentOmegaMinusInj;
	
	 	// cascades generated within acceptance (cut in pt + theta)
	 lHistThetaGenCascNat      = fHistThetaGenCascOmegaMinusNat;
         lHistThetaGenCascInj      = fHistThetaGenCascOmegaMinusInj;
	 l2dHistGenPtVsGenYFdbl = f2dHistGenPtVsGenYFdblOmegaMinus;
	 lHistThetaLambda       = fHistThetaLambdaOmegaMinus;
	 lHistThetaBach         = fHistThetaBachOmegaMinus;
	 lHistThetaBarDghter    = fHistThetaBarDghterOmegaMinus;
	 lHistThetaMesDghter    = fHistThetaMesDghterOmegaMinus;
	 lHistPtBach	        = fHistPtBachOmegaMinus;
	 lHistPtBarDghter       = fHistPtBarDghterOmegaMinus;
	 lHistPtMesDghter       = fHistPtMesDghterOmegaMinus;   
        break;
    
    case 4:  // Omega+
         lPdgCodeCasc       =  -3334;  //Omega+
         lPdgCodeBach       =    321;  //K+
         lPdgCodeLambda     =  -3122;  //AntiLambda0
         lPdgCodeDghtMesV0  =    211;  //Pi+
         lPdgCodeDghtBarV0  =  -2212;  //AntiProton 
	 
	 	// any Omega+
	 lHistEtaGenCasc        = fHistEtaGenCascOmegaPlus;
         l3dHistGenPtVsGenYvsCentNat = f3dHistGenPtVsGenYvsCentOmegaPlusNat;
         l3dHistGenPtVsGenYvsNtracksNat = f3dHistGenPtVsGenYvsNtracksOmegaPlusNat;
         l3dHistGenPtVsGenYvsCentInj = f3dHistGenPtVsGenYvsCentOmegaPlusInj;
         l3dHistGenPtVsGenYvsNtracksInj = f3dHistGenPtVsGenYvsNtracksOmegaPlusInj;
         l3dHistGenPtVsGenctauvsCentNat = f3dHistGenPtVsGenctauvsCentOmegaPlusNat;
         l3dHistGenPtVsGenctauvsCentInj = f3dHistGenPtVsGenctauvsCentOmegaPlusInj;

	 	
	 	// cascades generated within acceptance (cut in pt + theta)
	 lHistThetaGenCascNat      = fHistThetaGenCascOmegaPlusNat;
         lHistThetaGenCascInj      = fHistThetaGenCascOmegaPlusInj;
	 l2dHistGenPtVsGenYFdbl = f2dHistGenPtVsGenYFdblOmegaPlus;
	 lHistThetaLambda       = fHistThetaLambdaOmegaPlus;
	 lHistThetaBach         = fHistThetaBachOmegaPlus;
	 lHistThetaBarDghter    = fHistThetaBarDghterOmegaPlus;
	 lHistThetaMesDghter    = fHistThetaMesDghterOmegaPlus;
	 lHistPtBach	        = fHistPtBachOmegaPlus;
	 lHistPtBarDghter       = fHistPtBarDghterOmegaPlus;
	 lHistPtMesDghter       = fHistPtMesDghterOmegaPlus;  
        break;

  }// end switch cascade


  for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < nNumberOfMCPrimaries; iCurrentLabelStack++) {

    Bool_t kIsNaturalPart = kFALSE;
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
        //cout << "Xi- within loop " << iCurrentLabelStack << "/ " << nNumberOfMCPrimaries << endl;
		
	// -  Xi level ... _____________________________________________________________
	TParticle* xiMC = 0x0;
		   xiMC = lCurrentParticle;
	if (!xiMC) {
	  Printf("MC TParticle pointer to Cascade = 0x0 ! Skip ...");
	  continue;
	}

        kIsNaturalPart = lMCevent->IsFromBGEvent(iCurrentLabelStack); 
        if (!kIsNaturalPart) { 
          if (!(lCurrentParticle->GetFirstMother()<0)) kIsNaturalPart = kTRUE; // because there are primaries (ALICE definition) not produced in the collision   
        }
        partEnergy = xiMC->Energy();
        partPz = xiMC->Pz();
        partEta = xiMC->Eta();
        partPt = xiMC->Pt();
        partP = xiMC->P();
        partTheta = xiMC->Theta();
        partMass = xiMC->GetMass();
        partVx = xiMC->Vx();
        partVy = xiMC->Vy();
        partVz = xiMC->Vz();

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
      
      if (endOfHijingEvent==-1) {  // For injected MC: determine where HIJING event ends 
        if ((lCurrentParticleaod->GetStatus()==21)||
            ( (lCurrentParticleaod->GetPdgCode() == 443) &&
              (lCurrentParticleaod->GetMother() == -1)   &&
              (lCurrentParticleaod->GetDaughterLabel(0) ==  (iCurrentLabelStack+1)))) {
 
         
          endOfHijingEvent = iCurrentLabelStack; 

        }
      }  

      if (!lCurrentParticleaod->IsPhysicalPrimary()) continue;  
      if (!(lCurrentParticleaod->PdgCode() == lPdgCodeCasc)) continue;

       
      kIsNaturalPart = kTRUE;
      if ((iCurrentLabelStack>=endOfHijingEvent)&&(endOfHijingEvent!=-1)&&(lCurrentParticleaod->GetMother()<0)) kIsNaturalPart = kFALSE; 

      partEnergy = lCurrentParticleaod->E();
      partPz = lCurrentParticleaod->Pz();
      partEta = lCurrentParticleaod->Eta();
      partP = lCurrentParticleaod->P();
      partPt = lCurrentParticleaod->Pt();
      partTheta = lCurrentParticleaod->Theta();
      partMass = lCurrentParticleaod->M(); // FIXME not sure this works, seems not implemented
      partVx = lCurrentParticleaod->Xv();
      partVy = lCurrentParticleaod->Yv();
      partVz = lCurrentParticleaod->Zv();

      if (lCurrentParticleaod->GetDaughterLabel(0)>=0) {
        AliAODMCParticle *mcBach = (AliAODMCParticle*) arrayMC->At(lCurrentParticleaod->GetDaughterLabel(0));
        if (mcBach) {
          bacVx  = mcBach->Xv();
          bacVy  = mcBach->Yv();
          bacVz  = mcBach->Zv();
        }
      }

    }
      
    ncascperevtot++;	
	// Fill the first histos : = any generated Xi, not necessarily within the acceptance
    Double_t lRapXiMC = 0.5*TMath::Log((partEnergy + partPz) / (partEnergy - partPz +1.e-13));
        // Calculate proper time
    Double_t lctau = TMath::Sqrt((partVx-bacVx)*(partVx-bacVx)+(partVy-bacVy)*(partVy-bacVy)+(partVz-bacVz)*(partVz-bacVz));
    if (partP!=0.)    lctau = lctau*partMass/partP;
    else lctau = -1.;

    Double_t lRadToDeg = 180.0/TMath::Pi();

    // Fill the first histos : = any generated Xi, not necessarily within the acceptance		
    lHistEtaGenCasc	->Fill( partEta );	 

    if (kIsNaturalPart) {
      l3dHistGenPtVsGenYvsCentNat   ->Fill( partPt, lRapXiMC, lcentrality );
      l3dHistGenPtVsGenYvsNtracksNat->Fill( partPt, lRapXiMC, lPrimaryTrackMultiplicity );     
      if (TMath::Abs(lRapXiMC)<0.5) l3dHistGenPtVsGenctauvsCentNat->Fill( partPt, lctau, lcentrality );
      lHistThetaGenCascNat          ->Fill( lRadToDeg * partTheta );
    } else {
      l3dHistGenPtVsGenYvsCentInj   ->Fill( partPt, lRapXiMC, lcentrality );
      l3dHistGenPtVsGenYvsNtracksInj->Fill( partPt, lRapXiMC, lPrimaryTrackMultiplicity );
      if (TMath::Abs(lRapXiMC)<0.5) l3dHistGenPtVsGenctauvsCentInj->Fill( partPt, lctau, lcentrality );
      lHistThetaGenCascInj          ->Fill( lRadToDeg * partTheta );
    }

    // Check the emission of particle stays within the acceptance of the detector (cut in theta)
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
      if ( lDght0ofXi->GetPdgCode() == lPdgCodeLambda   &&  // Here !
           lDght1ofXi->GetPdgCode() == lPdgCodeBach ){      // Here !
				
        lLambda = lDght0ofXi;
        lBach   = lDght1ofXi;
      }// end if dghter 0 = Lambda and    dghter 1 = Pi-  
			
      // Xi - Case 2
        else if ( lDght0ofXi->GetPdgCode() == lPdgCodeBach  &&      // Here !
                  lDght1ofXi->GetPdgCode() == lPdgCodeLambda ){     // Here !
			
        lBach   = lDght0ofXi;
	lLambda = lDght1ofXi;
      }//  end if dghter 0 = Pi-  and   dghter 1 = Lambda
			
	// V0 otherwise - Case 3	
        else continue;
		
      // Check the emission of particle stays within the acceptance of the detector (cut in pt + theta)
      if (fApplyAccCut) { 
        if ( lLambda->Theta() < TMath::Pi()/4.0  ||    lLambda->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
        if( lBach->Theta() < TMath::Pi()/4.0    ||    lBach->Theta() > 3.0*TMath::Pi()/4.0 )   continue;
        if( lBach->Pt() < 0.150 ) continue; //FIXME : maybe tuned for Xi but not for K- from Omega ...
      } 		
		
      // -  V0 level ... _____________________________________________________________
      TParticle* lDghtBarV0 = 0;
      TParticle* lDghtMesV0 = 0;
	
      if( lLambda->GetNDaughters() != 2 )  continue;
      if( lLambda->GetDaughter(0) < 0 )    continue;
      if( lLambda->GetDaughter(1) < 0 )    continue;
		
		
      TParticle* lDght0ofLambda = lMCstack->Particle(  lLambda->GetDaughter(0) );
      TParticle* lDght1ofLambda = lMCstack->Particle(  lLambda->GetDaughter(1) );
			
      // V0 - Case 1
      if ( lDght0ofLambda->GetPdgCode() == lPdgCodeDghtBarV0 &&    // Here !
           lDght1ofLambda->GetPdgCode() == lPdgCodeDghtMesV0 ) {    // Here !
			
        lDghtBarV0 = lDght0ofLambda;
        lDghtMesV0 = lDght1ofLambda;
      }// end if dghter 0 = Proton  and   dghter 1 = Pi-  
			
      // V0 - Case 2
        else if ( lDght0ofLambda->GetPdgCode() == lPdgCodeDghtMesV0  &&     // Here !
                  lDght1ofLambda->GetPdgCode() == lPdgCodeDghtBarV0 ) {      // Here !
			
        lDghtMesV0 = lDght0ofLambda;
	lDghtBarV0 = lDght1ofLambda;
      }//  end if dghter 0 = Pi-  and   dghter 1 = proton
			
      // V0 otherwise - Case 3
        else continue;
	
			
      // Check the emission of particle stays within the acceptance of the detector
      if (fApplyAccCut) { 
        if( lDghtBarV0->Theta() < TMath::Pi()/4.0  ||  lDghtBarV0->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
        if( lDghtMesV0->Theta() < TMath::Pi()/4.0  ||  lDghtMesV0->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
	
        if( lDghtBarV0->Pt() < 0.250 ) continue;
        if( lDghtMesV0->Pt() < 0.150 ) continue;
      }

      lambdaTheta = lLambda->Theta();
      bacTheta    = lBach->Theta();
      dghtBarV0Theta = lDghtBarV0->Theta(); 			
      dghtMesV0Theta = lDghtMesV0->Theta();
      bacPt       = lBach->Pt();
      dghtBarV0Pt = lDghtBarV0->Pt();
      dghtMesV0Pt = lDghtMesV0->Pt();
	 	 	
    } else if ( fAnalysisType == "AOD") {

      AliAODMCParticle *xiMC = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack);
      if ( xiMC->GetNDaughters() != 2) continue;
      if ( xiMC->GetDaughterLabel(0) < 0 )  continue;
      if ( xiMC->GetDaughterLabel(1) < 0 )  continue;

      AliAODMCParticle* lDght0ofXi = (AliAODMCParticle*) arrayMC->At(  xiMC->GetDaughterLabel(0) );
      AliAODMCParticle* lDght1ofXi = (AliAODMCParticle*) arrayMC->At(  xiMC->GetDaughterLabel(1) );

      AliAODMCParticle* lLambda = 0;
      AliAODMCParticle* lBach   = 0;

      // Xi - Case 1
      if ( lDght0ofXi->PdgCode() == lPdgCodeLambda   &&  // Here !
           lDght1ofXi->PdgCode() == lPdgCodeBach ){      // Here !

        lLambda = lDght0ofXi;
        lBach   = lDght1ofXi;
      }// end if dghter 0 = Lambda and    dghter 1 = Pi-

      // Xi - Case 2
        else if ( lDght0ofXi->PdgCode() == lPdgCodeBach  &&      // Here !
                  lDght1ofXi->PdgCode() == lPdgCodeLambda ){     // Here !

        lBach   = lDght0ofXi;
        lLambda = lDght1ofXi;
      }//  end if dghter 0 = Pi-  and   dghter 1 = Lambda

        // V0 otherwise - Case 3
        else continue;

      // Check the emission of particle stays within the acceptance of the detector (cut in pt + theta)
      if (fApplyAccCut) {
        if ( lLambda->Theta() < TMath::Pi()/4.0  ||    lLambda->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
        if( lBach->Theta() < TMath::Pi()/4.0    ||    lBach->Theta() > 3.0*TMath::Pi()/4.0 )   continue;
        if( lBach->Pt() < 0.150 ) continue; //FIXME : maybe tuned for Xi but not for K- from Omega ...
      }

      // -  V0 level ... _____________________________________________________________
      AliAODMCParticle* lDghtBarV0 = 0;
      AliAODMCParticle* lDghtMesV0 = 0;

      if( lLambda->GetNDaughters() != 2 )  continue;
      if( lLambda->GetDaughterLabel(0) < 0 )    continue;
      if( lLambda->GetDaughterLabel(1) < 0 )    continue;


      AliAODMCParticle* lDght0ofLambda = (AliAODMCParticle*) arrayMC->At(  lLambda->GetDaughterLabel(0) );
      AliAODMCParticle* lDght1ofLambda = (AliAODMCParticle*) arrayMC->At(  lLambda->GetDaughterLabel(1) );

      // V0 - Case 1
      if ( lDght0ofLambda->PdgCode() == lPdgCodeDghtBarV0 &&    // Here !
           lDght1ofLambda->PdgCode() == lPdgCodeDghtMesV0 ) {    // Here !

        lDghtBarV0 = lDght0ofLambda;
        lDghtMesV0 = lDght1ofLambda;
      }// end if dghter 0 = Proton  and   dghter 1 = Pi-

      // V0 - Case 2
        else if ( lDght0ofLambda->PdgCode() == lPdgCodeDghtMesV0  &&     // Here !
                  lDght1ofLambda->PdgCode() == lPdgCodeDghtBarV0 ) {      // Here !

        lDghtMesV0 = lDght0ofLambda;
        lDghtBarV0 = lDght1ofLambda;
      }//  end if dghter 0 = Pi-  and   dghter 1 = proton

      // V0 otherwise - Case 3
        else continue;


      // Check the emission of particle stays within the acceptance of the detector
      if (fApplyAccCut) {
        if( lDghtBarV0->Theta() < TMath::Pi()/4.0  ||  lDghtBarV0->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
        if( lDghtMesV0->Theta() < TMath::Pi()/4.0  ||  lDghtMesV0->Theta() > 3.0*TMath::Pi()/4.0 ) continue;

        if( lDghtBarV0->Pt() < 0.250 ) continue;
        if( lDghtMesV0->Pt() < 0.150 ) continue;
      }


      lambdaTheta = lLambda->Theta();
      bacTheta    = lBach->Theta();
      dghtBarV0Theta = lDghtBarV0->Theta();
      dghtMesV0Theta = lDghtMesV0->Theta();
      bacPt       = lBach->Pt();
      dghtBarV0Pt = lDghtBarV0->Pt();
      dghtMesV0Pt = lDghtMesV0->Pt();

    }

    // - Just to know which file is currently open : locate the file containing Xi
                //cout << "Name of the file containing generated Xi :" <<  CurrentFileName()
                //       
    // - Filling histos for findable cascades... _________________________________________________________________

    // - Fill theta histos for Lambda and Bach
    lHistThetaLambda        ->Fill( lRadToDeg * lambdaTheta );
    lHistThetaBach          ->Fill( lRadToDeg * bacTheta );

    // - Fill theta histos for V0 daughters
    lHistThetaBarDghter     ->Fill( lRadToDeg * dghtBarV0Theta );
    lHistThetaMesDghter     ->Fill( lRadToDeg * dghtMesV0Theta );

    // - Fill pt histos.
    lHistPtBach             ->Fill( bacPt );
    lHistPtBarDghter        ->Fill( dghtBarV0Pt );
    lHistPtMesDghter        ->Fill( dghtMesV0Pt );
    //  if(iCascType == 1) Printf("Xi- current index = %n ", iCurrentLabelStack);

    l2dHistGenPtVsGenYFdbl  ->Fill( partPt, lRapXiMC );

    ncascperev++;			
		
	     
  }// This is the end of the loop on primaries
  
  if (iCascType == 1) {
    fHistnXiMinusPerEv->Fill(ncascperev);
    fHistnXiMinusPerEvTot->Fill(ncascperevtot);
//       Printf("N xi-tot = %n N xi-acc = %n\n", ncascperevtot, ncascperev);
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

   
// - Re-initialisation of the local TH1F pointers
lHistEtaGenCasc         = 0x0;
l3dHistGenPtVsGenYvsCentNat = 0x0;
l3dHistGenPtVsGenYvsNtracksNat = 0x0;
l3dHistGenPtVsGenYvsCentInj = 0x0;
l3dHistGenPtVsGenYvsNtracksInj = 0x0;
l3dHistGenPtVsGenctauvsCentNat = 0x0;
l3dHistGenPtVsGenctauvsCentInj = 0x0;

lHistThetaGenCascNat    = 0x0;
lHistThetaGenCascInj    = 0x0;
l2dHistGenPtVsGenYFdbl  = 0x0;
lHistThetaLambda        = 0x0;
lHistThetaBach          = 0x0;
lHistThetaBarDghter     = 0x0;
lHistThetaMesDghter     = 0x0;
lHistPtBach             = 0x0;
lHistPtBarDghter        = 0x0;
lHistPtMesDghter        = 0x0;	

} // end of loop over the different types of cascades (Xi-+, Omega-+)
 	
 
 
//__________________________________________________________________________	
// Part 3 - Loop over the reconstructed candidates
  
Int_t nAssoXiMinus = 0;
Int_t nAssoXiPlus = 0;
Int_t nAssoOmegaMinus = 0;
Int_t nAssoOmegaPlus = 0;


Int_t lPosTPCClusters   = 0;
Int_t lNegTPCClusters   = 0;
Int_t lBachTPCClusters  = 0;


// Double_t lChi2Xi         = -1. ;
Double_t lDcaXiDaughters            = -1. ;
Double_t lDcaBachToPrimVertexXi     = -1. ;
Double_t lXiCosineOfPointingAngle   = -1. ;
Double_t lPosXi[3]                  = { -1000.0, -1000.0, -1000.0 };
Double_t lXiRadius                  = -1000. ;

Double_t lInvMassLambdaAsCascDghter = 0.;
Double_t lDcaV0DaughtersXi          = -1.;
// Double_t lV0Chi2Xi               = -1. ;
Double_t lV0CosineOfPointingAngle = -1.;
Double_t lV0toXiCosineOfPointingAngle   = -1.;
Double_t lPosV0Xi[3]                = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
Double_t lV0RadiusXi                = -1000.;
Double_t lDcaV0ToPrimVertexXi       = -1.;
Double_t lDcaPosToPrimVertexXi      = -1.;
Double_t lDcaNegToPrimVertexXi      = -1.;

Double_t lChargeXi                  = -1.;

Double_t lV0mom = -1000.;
 
Double_t lmcPt = -1.;         
Double_t lmcRapCasc = -1.; 
Double_t lmcEta = -1000.;     
Double_t lmcTransvRadius = -1000.; 


Double_t lrecoPt = -100.;     
Double_t lrecoTransvRadius = -1000.; 

Double_t lDeltaPhiMcReco = -1.;


Double_t lBachTransvMom = 0.;
Double_t lpTrackTransvMom  = 0.;
Double_t lnTrackTransvMom  = 0.;

Double_t lmcPtPosV0Dghter = -100.;
Double_t lmcPtNegV0Dghter = -100.;
Double_t lrecoP = -100.;

Double_t   lmcPtBach = -100.;

Double_t cascadeMass = 0.;



for (Int_t iXi = 0; iXi < ncascades; iXi++) {// This is the begining of the Cascade loop


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

      // Combined PID
      // Reasonable guess for the priors for the cascade track sample (e-, mu, pi, K, p)
      Double_t lPriorsGuessXi[5]    = {0, 0, 2, 0, 1};
      Double_t lPriorsGuessOmega[5] = {0, 0, 1, 1, 1};

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

      Bool_t kIsNaturalPart = kFALSE;

      Float_t  etaBach = 0.;
      Float_t  etaPos  = 0.;
      Float_t  etaNeg  = 0.;


      if ( fAnalysisType == "ESD" ) {		
	AliESDcascade *xiESD = lESDevent->GetCascade(iXi);
	if (!xiESD) continue;
	
	// - Step II.1 : Connection to daughter tracks of the current cascade
	//-------------
			
		UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xiESD->GetPindex() );
		UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xiESD->GetNindex() );
		UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xiESD->GetBindex() );
		// abs value not needed ; the index should always be positive (!= label ...)
                
                
        // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
        if(lBachIdx == lIdxNegXi) {
                AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
        }
        if(lBachIdx == lIdxPosXi) {
                AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
        }
      
	AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
	AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
	AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx  );
	if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
		Printf("ERROR: Could not retrieve one of the 3 daughter tracks of the cascade ...");
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
        if(fkQualityCutnTPCcls){
                // 2 - Poor quality related to TPC clusters
                if(lPosTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Pos. track has less than n TPC clusters ... continue!"); continue; }
                if(lNegTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Neg. track has less than n TPC clusters ... continue!"); continue; }
                if(lBachTPCClusters < fMinnTPCcls) { AliWarning("Pb / Bach.   track has less than n TPC clusters ... continue!"); continue; }
        }

        etaPos  = pTrackXi->Eta();
        etaNeg  = nTrackXi->Eta();
        etaBach = bachTrackXi->Eta();
       
	
	// - Step II.2 : Info over reconstructed cascades
	//-------------	
	
	
	Double_t lV0quality = 0.;
	
	if( bachTrackXi->Charge() < 0 )	{
		lV0quality = 0.;
		xiESD->ChangeMassHypothesis(lV0quality , 3312); 	
			// Calculate the effective mass of the Xi- candidate. 
			// pdg code 3312 = Xi-
		lInvMassXiMinus = xiESD->GetEffMassXi();
		
		lV0quality = 0.;
		xiESD->ChangeMassHypothesis(lV0quality , 3334); 	
			// Calculate the effective mass of the Xi- candidate. 
			// pdg code 3334 = Omega-
		lInvMassOmegaMinus = xiESD->GetEffMassXi();
					
		lV0quality = 0.;
		xiESD->ChangeMassHypothesis(lV0quality , 3312); 	// Back to default hyp.
		
	}
	
	if( bachTrackXi->Charge() >  0 ){
		lV0quality = 0.;
		xiESD->ChangeMassHypothesis(lV0quality , -3312); 	
			// Calculate the effective mass of the Xi+ candidate. 
			// pdg code -3312 = Xi+
		lInvMassXiPlus = xiESD->GetEffMassXi();
		
		lV0quality = 0.;
		xiESD->ChangeMassHypothesis(lV0quality , -3334); 	
			// Calculate the effective mass of the Xi+ candidate. 
			// pdg code -3334  = Omega+
		lInvMassOmegaPlus = xiESD->GetEffMassXi();
		
		lV0quality = 0.;
		xiESD->ChangeMassHypothesis(lV0quality , -3312); 	// Back to "default" hyp.
	}
	

        //lChi2Xi                       = xiESD->GetChi2Xi();
        lDcaXiDaughters                 = xiESD->GetDcaXiDaughters();
        lDcaBachToPrimVertexXi          = TMath::Abs( bachTrackXi->GetD( lBestPrimaryVtxPos[0],
                                                                         lBestPrimaryVtxPos[1],
                                                                         lMagneticField  ) );
                                        // NOTE : AliExternalTrackParam::GetD returns an algebraic value
        lXiCosineOfPointingAngle        = xiESD->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                                  lBestPrimaryVtxPos[1],
                                                                                  lBestPrimaryVtxPos[2] );
                                        // Take care : the best available vertex should be used (like in AliCascadeVertexer)
        xiESD->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] );  
        lInvMassLambdaAsCascDghter      = xiESD->GetEffMass();
                                        // This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
        lDcaV0DaughtersXi               = xiESD->GetDcaV0Daughters();
        // lV0Chi2Xi                    = xiESD->GetChi2V0();
        lV0toXiCosineOfPointingAngle    = xiESD->GetV0CosineOfPointingAngle( lPosXi[0],    
                                                                             lPosXi[1],
                                                                             lPosXi[2] );

        lV0CosineOfPointingAngle      = xiESD->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                           lBestPrimaryVtxPos[1],
                                                                           lBestPrimaryVtxPos[2]); 

        //if (lV0toXiCosineOfPointingAngle==1.) cout << "Cosine V0 PA wrt Xi pos ==1!" <<endl;
        //else if (lV0toXiCosineOfPointingAngle>1.) cout <<"Cosine V0 PA wrt Xi pos >1!"<<endl;
        xiESD->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] );

        lDcaV0ToPrimVertexXi            = xiESD->GetD( lBestPrimaryVtxPos[0],
                                                lBestPrimaryVtxPos[1],
                                                lBestPrimaryVtxPos[2] );

        lDcaPosToPrimVertexXi           = TMath::Abs( pTrackXi   ->GetD( lBestPrimaryVtxPos[0],
                                                                         lBestPrimaryVtxPos[1],
                                                                         lMagneticField  )     );

        lDcaNegToPrimVertexXi           = TMath::Abs( nTrackXi   ->GetD( lBestPrimaryVtxPos[0],
                                                                         lBestPrimaryVtxPos[1],
                                                                         lMagneticField  )     );

	lChargeXi = xiESD->Charge();
	
	// - Step II.3 : PID info   

	//-------------
	
	
	// 3.1 - PID Information


	// Combined VO-positive-daughter PID
	AliPID pPidXi;		pPidXi.SetPriors(    lPriorsGuessXi    );
	AliPID pPidOmega;	pPidOmega.SetPriors( lPriorsGuessOmega );
		
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
	AliPID nPidXi;		nPidXi.SetPriors(    lPriorsGuessXi    );
	AliPID nPidOmega;	nPidOmega.SetPriors( lPriorsGuessOmega );
		
	if( nTrackXi->IsOn(AliESDtrack::kESDpid) ) {  // Combined PID exists
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
	AliPID bachPidXi;	bachPidXi.SetPriors(    lPriorsGuessXi    );
	AliPID bachPidOmega;	bachPidOmega.SetPriors( lPriorsGuessOmega );
	
        
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
	
	
	// 3.1.B - TPC PID : 4-sigma bands on Bethe-Bloch curve

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
                
        }
*/        
		
	// 3.2 - PID proba Vs Pt(Bach)
	Int_t      lblBachForPID  = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
	TParticle* mcBachForPID   = lMCstack->Particle( lblBachForPID );
	lmcPtBach      = mcBachForPID->Pt();
	
	// 3.3 - MC perfect PID
	
	if( mcBachForPID->GetPdgCode() == -211) lIsBachelorMCPiMinus = kTRUE;
	if( mcBachForPID->GetPdgCode() ==  211) lIsBachelorMCPiPlus  = kTRUE;
	if( mcBachForPID->GetPdgCode() == -321) lIsBachelorMCKMinus  = kTRUE;
	if( mcBachForPID->GetPdgCode() ==  321) lIsBachelorMCKPlus   = kTRUE;
	
	
	
	// - Step II.4 : MC association (care : lots of "continue;" below this line)
	//-------------	
	
	if(fDebug > 5)
		cout 	<< "MC EventNumber : " << lMCevent->Header()->GetEvent() 
			<< " / MC event Number in Run : " << lMCevent->Header()->GetEventNrInRun() << endl;
	
	// - Step 4.1 : level of the V0 daughters

	Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );  
		// Abs value = needed ! question of quality track association ... (negative label when at least one cluster in the track is from a different particle)
	Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );  		
	TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
	TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
	

	// - Step 4.2 : level of the Xi daughters
		
	Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ; 
	Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();
	
		if( lblMotherPosV0Dghter != lblMotherNegV0Dghter) continue; // same mother
		if( lblMotherPosV0Dghter < 0 ) continue; // this particle is primary, no mother   
		if( lblMotherNegV0Dghter < 0 ) continue;
					

		// mothers = Lambda candidate ... a priori
	
	TParticle* mcMotherPosV0Dghter = lMCstack->Particle( lblMotherPosV0Dghter );
	TParticle* mcMotherNegV0Dghter = lMCstack->Particle( lblMotherNegV0Dghter );  // MN: redundant?? already checked that labels are the same...-->same part from stack

	Int_t      lblBach  = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
	TParticle* mcBach   = lMCstack->Particle( lblBach );	
				

	// - Step 4.3 : level of Xi candidate
	
	Int_t lblGdMotherPosV0Dghter =   mcMotherPosV0Dghter->GetFirstMother() ;
	Int_t lblGdMotherNegV0Dghter =   mcMotherNegV0Dghter->GetFirstMother() ;
			
		if( lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter ) continue;
		if( lblGdMotherPosV0Dghter < 0 ) continue; // primary lambda ...   
		if( lblGdMotherNegV0Dghter < 0 ) continue; // primary lambda ...   			
		
		// Gd mothers = Xi candidate ... a priori
	
	TParticle* mcGdMotherPosV0Dghter = lMCstack->Particle( lblGdMotherPosV0Dghter );
	TParticle* mcGdMotherNegV0Dghter = lMCstack->Particle( lblGdMotherNegV0Dghter );
					
	Int_t lblMotherBach = (Int_t) TMath::Abs( mcBach->GetFirstMother()  );  
	
	if( lblMotherBach != lblGdMotherPosV0Dghter ) continue; //same mother for bach and V0 daughters
	
	TParticle* mcMotherBach = lMCstack->Particle( lblMotherBach );
	
        // Check if cascade is primary

        if (!(lMCstack->IsPhysicalPrimary(lblMotherBach))) continue;  

        kIsNaturalPart = lMCevent->IsFromBGEvent(lblMotherBach); // if kTRUE the particle is a natural one
        if (!kIsNaturalPart) {
          if (!(mcMotherBach->GetFirstMother()<0)) kIsNaturalPart = kTRUE;
        }

	// - Step 4.4 : Manage boolean for association
    	
	if( mcMotherBach 		->GetPdgCode() ==   3312 &&
	    mcGdMotherPosV0Dghter	->GetPdgCode() ==   3312 &&
	    mcGdMotherNegV0Dghter	->GetPdgCode() ==   3312)    {	lAssoXiMinus = kTRUE;
                                                                        cascadeMass = 1.321;
                                                                      //  Printf("Xi- asso current index = %n ", lblGdMotherPosV0Dghter);  
                                                                        nAssoXiMinus++; }
	
	else if( mcMotherBach 		->GetPdgCode() ==  -3312 &&
	    mcGdMotherPosV0Dghter	->GetPdgCode() ==  -3312 &&
	    mcGdMotherNegV0Dghter	->GetPdgCode() ==  -3312)    {	lAssoXiPlus = kTRUE;
                                                                        cascadeMass = 1.321;
                                                                        nAssoXiPlus++; }
	
	else if( mcMotherBach 		->GetPdgCode() ==   3334 &&
	    mcGdMotherPosV0Dghter	->GetPdgCode() ==   3334 &&
	    mcGdMotherNegV0Dghter	->GetPdgCode() ==   3334)    {	lAssoOmegaMinus = kTRUE;
                                                                        cascadeMass = 1.672;
                                                                        nAssoOmegaMinus++; }
		
	else if( mcMotherBach 		->GetPdgCode() ==  -3334 &&
	    mcGdMotherPosV0Dghter	->GetPdgCode() ==  -3334 &&
	    mcGdMotherNegV0Dghter	->GetPdgCode() ==  -3334)    { 	lAssoOmegaPlus = kTRUE;
                                                                        cascadeMass = 1.672;
                                                                        nAssoOmegaPlus++; }
	
	
	
	// If a proper association  exists ...
		
	if(fDebug > 4){
		cout << "XiMinus    = " << lAssoXiMinus    << endl;
		cout << "XiPlus     = " << lAssoXiPlus     << endl;
		cout << "OmegaMinus = " << lAssoOmegaMinus << endl;
		cout << "OmegaPlus  = " << lAssoOmegaPlus  << endl 
		     << "----" 		<< endl;	
	}


	if(fDebug > 5){
		cout << endl;
		cout << "- V0 daughters - " << endl;
		cout << "     + V0 Pos. / Label : " << lblPosV0Dghter 
		<< " - Pdg Code : " << mcPosV0Dghter->GetTitle() << endl;
		cout << "     - V0 Neg. / Label : " << lblNegV0Dghter 
		<< " - Pdg Code : " << mcNegV0Dghter->GetTitle() << endl;
		
		cout << "- Xi daughters - " << endl;
		cout << "     + V0 Pos. mother / Label : " << lblMotherPosV0Dghter 
		<< " - Pdg Code : " << mcMotherPosV0Dghter->GetTitle() << endl;
		cout << "     - V0 Neg. mother / Label : " << lblMotherNegV0Dghter 
		<< " - Pdg Code : " << mcMotherNegV0Dghter->GetTitle() << endl;
		
		cout << "     --  Bach. / Label :" << lblBach 
		<< " -  Pdg Code : " << mcBach->GetTitle() << endl;
		
		cout << "- Xi candidate -" << endl;
		cout << "    +  V0 Pos. Gd Mother / Label : " << lblGdMotherPosV0Dghter 
		<< " - Pdg Code : " << mcGdMotherPosV0Dghter->GetTitle() << endl;
		cout << "    -  V0 Neg. Gd Mother / Label : "  << lblGdMotherNegV0Dghter 
		<< " - Pdg Code : "<< mcGdMotherNegV0Dghter->GetTitle() << endl;
		
		cout << "    --  Mother Bach. / Label : " << lblMotherBach 
		<< " - Pdg Code    : " << mcMotherBach->GetTitle() << endl;
		cout << endl;
	}

	
        lmcPt             = mcMotherBach->Pt();
        lmcRapCasc        = 0.5*TMath::Log( (mcMotherBach->Energy() + mcMotherBach->Pz()) / 
                                                     (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) );
        lmcEta            = mcMotherBach->Eta();
        lmcTransvRadius   = mcBach->R(); // to get the decay point of Xi, = the production vertex of Bachelor ...
        
        TVector3 lmcTVect3Mom( mcMotherBach->Px(), mcMotherBach->Py(), mcMotherBach->Pz() );

        lrecoPt           = xiESD->Pt();
        lrecoTransvRadius = TMath::Sqrt( xiESD->Xv() * xiESD->Xv() + xiESD->Yv() * xiESD->Yv() );
        
        TVector3 lrecoTVect3Mom( xiESD->Px(), xiESD->Py(), xiESD->Pz() );
        lDeltaPhiMcReco   = lmcTVect3Mom.DeltaPhi( lrecoTVect3Mom ) * 180.0/TMath::Pi();

        lmcPtPosV0Dghter = mcPosV0Dghter->Pt() ;
        lmcPtNegV0Dghter = mcNegV0Dghter->Pt();
        lrecoP           = xiESD->P();

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

        const AliAODcascade *xiAOD = lAODevent->GetCascade(iXi);

        if (!xiAOD) continue;

        // - Step II.1 : Connection to daughter tracks of the current cascade
        //-------------

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

                // abs value not needed ; the index should always be positive (!= label ...)


        // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
        if(lBachIdx == lIdxNegXi) {
                AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
        }
        if(lBachIdx == lIdxPosXi) {
                AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
        }

        lPosTPCClusters   = pTrackXi->GetTPCNcls();
        lNegTPCClusters   = nTrackXi->GetTPCNcls();
        lBachTPCClusters  = bachTrackXi->GetTPCNcls();


        // FIXME : rejection of a poor quality tracks
        if (fkQualityCutTPCrefit) {
                // 1 - Poor quality related to TPCrefit
          if (!(pTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
          if (!(nTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
          if (!(bachTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }

        }
        if (fkQualityCutnTPCcls) {
                // 2 - Poor quality related to TPC clusters
                if(lPosTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Pos. track has less than n TPC clusters ... continue!"); continue; }
                if(lNegTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Neg. track has less than n TPC clusters ... continue!"); continue; }
                if(lBachTPCClusters < fMinnTPCcls) { AliWarning("Pb / Bach.   track has less than n TPC clusters ... continue!"); continue; }
        }

        etaPos  = pTrackXi->Eta();
        etaNeg  = nTrackXi->Eta();
        etaBach = bachTrackXi->Eta();

        // - Step II.2 : Info over reconstructed cascades
        //-------------


        if( bachTrackXi->Charge() < 0 ) {
                lInvMassXiMinus = xiAOD->MassXi();

                lInvMassOmegaMinus = xiAOD->MassOmega();
        }

        if( bachTrackXi->Charge() >  0 ){
                lInvMassXiPlus = xiAOD->MassXi();

                lInvMassOmegaPlus = xiAOD->MassOmega();

        }


        //lChi2Xi                       = xiAOD->Chi2Xi();
        lDcaXiDaughters                 = xiAOD->DcaXiDaughters();
        lDcaBachToPrimVertexXi          = xiAOD->DcaBachToPrimVertex();
        lXiCosineOfPointingAngle        = xiAOD->CosPointingAngleXi( lBestPrimaryVtxPos[0],
                                                                     lBestPrimaryVtxPos[1],
                                                                     lBestPrimaryVtxPos[2] );

        lPosXi[0] = xiAOD->DecayVertexXiX();
        lPosXi[1] = xiAOD->DecayVertexXiY();
        lPosXi[2] = xiAOD->DecayVertexXiZ();

                                        // This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
        lDcaV0DaughtersXi               = xiAOD->DcaV0Daughters();
        // lV0Chi2Xi                    = xiAOD->GetChi2V0();
        lV0toXiCosineOfPointingAngle    = xiAOD->CosPointingAngle( lPosXi );

        lV0CosineOfPointingAngle        = xiAOD->CosPointingAngle( lBestPrimaryVtxPos );

        lPosV0Xi[0] = xiAOD->DecayVertexV0X();
        lPosV0Xi[1] = xiAOD->DecayVertexV0Y();
        lPosV0Xi[2] = xiAOD->DecayVertexV0Z();


        lDcaV0ToPrimVertexXi            = xiAOD->DcaV0ToPrimVertex();

        lDcaPosToPrimVertexXi           = xiAOD->DcaPosToPrimVertex();

        lDcaNegToPrimVertexXi           = xiAOD->DcaNegToPrimVertex();

        lChargeXi = xiAOD->ChargeXi();

        // - Step II.3 : PID info

        //-------------


        // 3.1 - PID Information


        // Combined VO-positive-daughter PID

        // Combined bachelor PID
/*        AliPID bachPidXi;       bachPidXi.SetPriors(    lPriorsGuessXi    );
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

        // 3.1.B - TPC PID : 4-sigma bands on Bethe-Bloch curve


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

        }
*/

        // 3.2 - PID proba Vs Pt(Bach)
        Int_t      lblBachForPID  = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
        AliAODMCParticle* mcBachForPID   = (AliAODMCParticle*) arrayMC->At( lblBachForPID );
        lmcPtBach      = mcBachForPID->Pt();

        // 3.3 - MC perfect PID

        if( mcBachForPID->PdgCode() == -211) lIsBachelorMCPiMinus = kTRUE;
        if( mcBachForPID->PdgCode() ==  211) lIsBachelorMCPiPlus  = kTRUE;
        if( mcBachForPID->PdgCode() == -321) lIsBachelorMCKMinus  = kTRUE;
        if( mcBachForPID->PdgCode() ==  321) lIsBachelorMCKPlus   = kTRUE;



        // - Step II.4 : MC association (care : lots of "continue;" below this line)
        //-------------

        if(fDebug > 5)
                cout    << "MC EventNumber : " << lMCevent->Header()->GetEvent()
                        << " / MC event Number in Run : " << lMCevent->Header()->GetEventNrInRun() << endl;

        // - Step 4.1 : level of the V0 daughters

        Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );
                // Abs value = needed ! question of quality track association ... (negative label when at least one cluster in the track is from a different particle)
        Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );
        AliAODMCParticle* mcPosV0Dghter = (AliAODMCParticle*) arrayMC->At( lblPosV0Dghter );
        AliAODMCParticle* mcNegV0Dghter = (AliAODMCParticle*) arrayMC->At( lblNegV0Dghter );


        // - Step 4.2 : level of the Xi daughters

        Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetMother();   
        Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetMother();

                if( lblMotherPosV0Dghter != lblMotherNegV0Dghter) continue; // same mother
                if( lblMotherPosV0Dghter < 0 ) continue; // this particle is primary, no mother
                if( lblMotherNegV0Dghter < 0 ) continue;


                // mothers = Lambda candidate ... a priori

        AliAODMCParticle* mcMotherPosV0Dghter = (AliAODMCParticle*) arrayMC->At( lblMotherPosV0Dghter );
        AliAODMCParticle* mcMotherNegV0Dghter = (AliAODMCParticle*) arrayMC->At( lblMotherNegV0Dghter );  // MN: redundant?? already checked that labels are the same...-->same part from stack

        Int_t      lblBach  = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
        AliAODMCParticle* mcBach   = (AliAODMCParticle*) arrayMC->At( lblBach );


        // - Step 4.3 : level of Xi candidate

        Int_t lblGdMotherPosV0Dghter =   mcMotherPosV0Dghter->GetMother() ;
        Int_t lblGdMotherNegV0Dghter =   mcMotherNegV0Dghter->GetMother() ;

                if( lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter ) continue;
                if( lblGdMotherPosV0Dghter < 0 ) continue; // primary lambda ...
                if( lblGdMotherNegV0Dghter < 0 ) continue; // primary lambda ...

                // Gd mothers = Xi candidate ... a priori

        AliAODMCParticle* mcGdMotherPosV0Dghter = (AliAODMCParticle*) arrayMC->At( lblGdMotherPosV0Dghter );
        AliAODMCParticle* mcGdMotherNegV0Dghter = (AliAODMCParticle*) arrayMC->At( lblGdMotherNegV0Dghter );

        Int_t lblMotherBach = (Int_t) TMath::Abs( mcBach->GetMother() );

        if( lblMotherBach != lblGdMotherPosV0Dghter ) continue; //same mother for bach and V0 daughters

        AliAODMCParticle* mcMotherBach = (AliAODMCParticle*) arrayMC->At( lblMotherBach );

                // Check if cascade is primary

        if (!(mcMotherBach->IsPhysicalPrimary())) continue;


        kIsNaturalPart = kTRUE;

        if ((lblMotherBach>=endOfHijingEvent)&&(endOfHijingEvent!=-1)&&(mcMotherBach->GetMother()<0)) kIsNaturalPart = kFALSE;

        // - Step 4.4 : Manage boolean for association

        if( mcMotherBach                ->GetPdgCode() ==   3312 &&
            mcGdMotherPosV0Dghter       ->GetPdgCode() ==   3312 &&
            mcGdMotherNegV0Dghter       ->GetPdgCode() ==   3312)    {  lAssoXiMinus = kTRUE;
                                                                        lInvMassLambdaAsCascDghter = xiAOD->MassLambda();
                                                                        cascadeMass = 1.321;
                                                                      //  Printf("Xi- asso current index = %n ", lblGdMotherPosV0Dghter);
                                                                        nAssoXiMinus++; }

        else if( mcMotherBach           ->GetPdgCode() ==  -3312 &&
            mcGdMotherPosV0Dghter       ->GetPdgCode() ==  -3312 &&
            mcGdMotherNegV0Dghter       ->GetPdgCode() ==  -3312)    {  lAssoXiPlus = kTRUE;
                                                                        lInvMassLambdaAsCascDghter = xiAOD->MassAntiLambda();
                                                                        cascadeMass = 1.321;
                                                                        nAssoXiPlus++; }

        else if( mcMotherBach           ->GetPdgCode() ==   3334 &&
            mcGdMotherPosV0Dghter       ->GetPdgCode() ==   3334 &&
            mcGdMotherNegV0Dghter       ->GetPdgCode() ==   3334)    {  lAssoOmegaMinus = kTRUE;
                                                                        lInvMassLambdaAsCascDghter = xiAOD->MassLambda();
                                                                        cascadeMass = 1.672;
                                                                        nAssoOmegaMinus++; }

        else if( mcMotherBach           ->GetPdgCode() ==  -3334 &&
            mcGdMotherPosV0Dghter       ->GetPdgCode() ==  -3334 &&
            mcGdMotherNegV0Dghter       ->GetPdgCode() ==  -3334)    {  lAssoOmegaPlus = kTRUE;
                                                                        lInvMassLambdaAsCascDghter = xiAOD->MassAntiLambda();
                                                                        cascadeMass = 1.672;
                                                                        nAssoOmegaPlus++; }
        //-------------

        lmcPt             = mcMotherBach->Pt();
        lmcRapCasc        = 0.5*TMath::Log( (mcMotherBach->E() + mcMotherBach->Pz()) /
                                                     (mcMotherBach->E() - mcMotherBach->Pz() +1.e-13) );
        lmcEta            = mcMotherBach->Eta();
        Float_t decayCascX = mcBach->Xv();
        Float_t decayCascY = mcBach->Yv();
        lmcTransvRadius   = TMath::Sqrt(decayCascX*decayCascX+decayCascY*decayCascY); // decay point of Xi, = the production vertex of Bachelor ...

        TVector3 lmcTVect3Mom( mcMotherBach->Px(), mcMotherBach->Py(), mcMotherBach->Pz() );

        Double_t xiMomX = xiAOD->MomXiX();
        Double_t xiMomY = xiAOD->MomXiY();
        Double_t xiMomZ = xiAOD->MomXiZ();

        lrecoPt           = TMath::Sqrt( xiMomX*xiMomX   + xiMomY*xiMomY ); 
        lrecoTransvRadius = TMath::Sqrt( xiAOD->DecayVertexXiX() * xiAOD->DecayVertexXiX() + xiAOD->DecayVertexXiY() * xiAOD->DecayVertexXiY() );

        TVector3 lrecoTVect3Mom( xiMomX, xiMomY, xiMomZ );
        lDeltaPhiMcReco   = lmcTVect3Mom.DeltaPhi( lrecoTVect3Mom ) * 180.0/TMath::Pi();

        lmcPtPosV0Dghter = mcPosV0Dghter->Pt() ;
        lmcPtNegV0Dghter = mcNegV0Dghter->Pt();
        lrecoP           = TMath::Sqrt( xiMomX*xiMomX   + xiMomY*xiMomY   + xiMomZ*xiMomZ );;

        Double_t lV0momX = xiAOD->MomV0X();
        Double_t lV0momY = xiAOD->MomV0Y();
        Double_t lV0momZ = xiAOD->MomV0Z();
        lV0mom = TMath::Sqrt(TMath::Power(lV0momX,2)+TMath::Power(lV0momY,2)+TMath::Power(lV0momZ,2));


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

      lXiRadius                       = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
      lV0RadiusXi                     = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] ); 


      // Cut on pt of the three daughter tracks
      if (lBachTransvMom<fMinPtCutOnDaughterTracks) continue;
      if (lpTrackTransvMom<fMinPtCutOnDaughterTracks) continue;
      if (lnTrackTransvMom<fMinPtCutOnDaughterTracks) continue;

      // Cut on pseudorapidity of the three daughter tracks
      if (TMath::Abs(etaBach)>fEtaCutOnDaughterTracks) continue;
      if (TMath::Abs(etaPos)>fEtaCutOnDaughterTracks) continue;
      if (TMath::Abs(etaNeg)>fEtaCutOnDaughterTracks) continue;


      // Extra-selection for cascade candidates
                // Towards optimisation of AA selection
      if (fkExtraSelections) {

        // if(lChi2Xi > 2000) continue; // in AliCascadeVertexer
        // if(lV0Chi2Xi > 2000) continue; // in AliV0vertexer

        if (lDcaXiDaughters > 0.3) continue; // in AliCascadeVertexer
        if (lXiCosineOfPointingAngle < 0.999 ) continue; // in AliCascadeVertexer
        if (lDcaV0ToPrimVertexXi < 0.05) continue; // in AliCascadeVertexer
        if (lDcaBachToPrimVertexXi < 0.03) continue; // in AliCascadeVertexer
////      if (TMath::Abs(lInvMassLambdaAsCascDghter-1.11568) > 0.006 ) continue;  // in AliCascadeVertexer

        if (lDcaV0DaughtersXi > 1.) continue; // in AliV0vertexer
        if (lV0toXiCosineOfPointingAngle < 0.998) continue; // in AliV0vertexer
        if (lDcaPosToPrimVertexXi < 0.1) continue; // in AliV0vertexer
        if (lDcaNegToPrimVertexXi < 0.1) continue; // in AliV0vertexer


          if(lXiRadius < .9) continue; // in AliCascadeVertexer
//        if(lXiRadius > 100) continue; // in AliCascadeVertexer
          if(lV0RadiusXi < 0.9) continue; // in AliV0vertexer
//        if(lV0RadiusXi > 100) continue; // in AliV0vertexer

      }


      // Combined PID TH1s
      if( lChargeXi < 0 && lIsBachelorPion )    fHistMassWithCombPIDXiMinus     ->Fill( lInvMassXiMinus    );
      if( lChargeXi > 0 && lIsBachelorPion )    fHistMassWithCombPIDXiPlus      ->Fill( lInvMassXiPlus     );
      if( lChargeXi < 0 && lIsBachelorKaon )    fHistMassWithCombPIDOmegaMinus  ->Fill( lInvMassOmegaMinus );
      if( lChargeXi > 0 && lIsBachelorKaon )    fHistMassWithCombPIDOmegaPlus   ->Fill( lInvMassOmegaPlus  );

 
      if( lChargeXi < 0 )     fHistMassXiMinus        ->Fill( lInvMassXiMinus );
      if( lChargeXi > 0 )     fHistMassXiPlus         ->Fill( lInvMassXiPlus );
      if( lChargeXi < 0 )     fHistMassOmegaMinus     ->Fill( lInvMassOmegaMinus );
      if( lChargeXi > 0 )     fHistMassOmegaPlus      ->Fill( lInvMassOmegaPlus );

      if(lIsBachelorPion)   f2dHistPIDprobaPionVsMCPtBach->Fill( lmcPtBach, ppionBach );
      if(lIsBachelorKaon)   f2dHistPIDprobaKaonVsMCPtBach->Fill( lmcPtBach, pkaonBach );

      if( lChargeXi < 0 && lIsBachelorMCPiMinus )    fHistMassWithMcPIDXiMinus     ->Fill( lInvMassXiMinus );
      if( lChargeXi > 0 && lIsBachelorMCPiPlus  )    fHistMassWithMcPIDXiPlus      ->Fill( lInvMassXiPlus );
      if( lChargeXi < 0 && lIsBachelorMCKMinus  )    fHistMassWithMcPIDOmegaMinus  ->Fill( lInvMassOmegaMinus );
      if( lChargeXi > 0 && lIsBachelorMCKPlus   )    fHistMassWithMcPIDOmegaPlus   ->Fill( lInvMassOmegaPlus );



      if(!lAssoXiMinus && !lAssoXiPlus && !lAssoOmegaMinus && !lAssoOmegaPlus) continue; // no association, skip the rest of the code
        
      // - Fill pt histos. xim only!!
      if (lAssoXiMinus) {
        fHistPtRecBachXiMinus->Fill(lBachTransvMom);
        fHistPtRecMesDghterXiMinus->Fill(lpTrackTransvMom);
        fHistPtRecBarDghterXiMinus->Fill(lnTrackTransvMom);
      }
 
        // Calculate proper time for cascade (reconstructed)
         
        Double_t lctau =  TMath::Sqrt(TMath::Power((lPosXi[0]-lBestPrimaryVtxPos[0]),2)+TMath::Power((lPosXi[1]-lBestPrimaryVtxPos[1]),2)+TMath::Power((lPosXi[2]-lBestPrimaryVtxPos[2]),2));
        if (lrecoP!=0)         lctau = lctau*cascadeMass/lrecoP;   
        else lctau = -1.;


        // Calculate proper time for Lambda (reconstructed)
        Float_t lambdaMass = 1.115683; // PDG mass 
        Float_t distV0Xi =  TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2)+TMath::Power((lPosV0Xi[2]-lPosXi[2]),2)); 
        Float_t lctauV0 = -1.;
        if (lV0mom!=0) lctauV0 = distV0Xi*lambdaMass/lV0mom; 

        Float_t distTV0Xi =  TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2));
       
 
	// - Histos for the cascade candidates associated with MC
	
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

        fHistV0toXiCosineOfPointingAngle->Fill(lV0toXiCosineOfPointingAngle);
       
        // - Step 6 : Containers = Cascade cuts + PID
	//-------------	
        // 6.3 - Filling the AliCFContainer (optimisation of topological selections + systematics)
        Double_t lContainerCutVars[22] = {0.0};

        lContainerCutVars[0]  = lDcaXiDaughters;
        lContainerCutVars[1]  = lDcaBachToPrimVertexXi;
        lContainerCutVars[2]  = lXiCosineOfPointingAngle;
        lContainerCutVars[3]  = lXiRadius;
        lContainerCutVars[4]  = lInvMassLambdaAsCascDghter;
        lContainerCutVars[5]  = lDcaV0DaughtersXi;
        lContainerCutVars[6]  = lV0CosineOfPointingAngle;
        lContainerCutVars[7]  = lV0RadiusXi;
        lContainerCutVars[8]  = lDcaV0ToPrimVertexXi;	
        lContainerCutVars[9]  = lDcaPosToPrimVertexXi;
        lContainerCutVars[10] = lDcaNegToPrimVertexXi;

        lContainerCutVars[13] = lmcPt;

        if (kIsNaturalPart) lContainerCutVars[16] = 0.;
        else lContainerCutVars[16] = 1.;
        lContainerCutVars[17] = lcentrality;
        lContainerCutVars[18] = lPrimaryTrackMultiplicity;       
        lContainerCutVars[19] = lctau;
        lContainerCutVars[20] = lctauV0;
        lContainerCutVars[21] = distTV0Xi;

        // All cases should be covered below
        if( lChargeXi < 0 && lAssoXiMinus    ) {
                lContainerCutVars[11] = lInvMassXiMinus;
                lContainerCutVars[12] = lInvMassOmegaMinus;//1.63;
                lContainerCutVars[14] = lmcRapCasc;
                lContainerCutVars[15] = -1.;
                if ( lIsBachelorPionForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC ) {   
                  fCFContAsCascadeCuts->Fill(lContainerCutVars,0); // for Xi-
                  fHistEtaBachXiM->Fill(etaBach);
                  fHistEtaPosXiM->Fill(etaPos);
                  fHistEtaNegXiM->Fill(etaNeg);
                }
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
        
        
	// 6.4 - Filling the AliCFContainers related to PID

	Double_t lContainerPIDVars[4] = {0.0};

	
	// Xi Minus		
	if( lChargeXi < 0 && lAssoXiMinus ) {
		lContainerPIDVars[0] = lmcPt              ;
		lContainerPIDVars[1] = lInvMassXiMinus    ;
		lContainerPIDVars[2] = lmcRapCasc         ;
		lContainerPIDVars[3] = lcentrality;
			
		// No PID
			fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorPionForTPC  )
			fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
		if( lIsBachelorPionForTPC && 
		    lIsPosProtonForTPC     )
			fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorPionForTPC && 
		    lIsPosProtonForTPC    && 
		    lIsNegPionForTPC       )
			fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorPion        )
			fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorPion       && 
		    lIsPosInXiProton    )
			fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if( lIsBachelorPion     && 
		    lIsPosInXiProton && 
		    lIsNegInXiPion    )
		 	fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Xi Plus		
	if( lChargeXi > 0 && lAssoXiPlus ) {
		lContainerPIDVars[0] = lmcPt         ;
		lContainerPIDVars[1] = lInvMassXiPlus;
		lContainerPIDVars[2] = lmcRapCasc    ;
		lContainerPIDVars[3] = lcentrality   ;
			
		// No PID
			fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorPionForTPC  )
			fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
		if( lIsBachelorPionForTPC && 
		    lIsNegProtonForTPC     )
			fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorPionForTPC && 
		    lIsNegProtonForTPC    && 
		    lIsPosPionForTPC       )
			fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorPion        )
			fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorPion       && 
		    lIsNegInXiProton    )
			fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if( lIsBachelorPion     && 
		    lIsNegInXiProton && 
		    lIsPosInXiPion    )
		 	fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Omega Minus		
	if( lChargeXi < 0 && lAssoOmegaMinus ) {
		lContainerPIDVars[0] = lmcPt              ;
		lContainerPIDVars[1] = lInvMassOmegaMinus ;
		lContainerPIDVars[2] = lmcRapCasc         ;
		lContainerPIDVars[3] = lcentrality;
			
		// No PID
			fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorKaonForTPC  )
			fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
		if( lIsBachelorKaonForTPC && 
		    lIsPosProtonForTPC     )
			fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorKaonForTPC && 
		    lIsPosProtonForTPC    && 
		    lIsNegPionForTPC       )
			fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorKaon        )
			fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorKaon       && 
		    lIsPosInOmegaProton    )
			fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if(lIsBachelorKaon     && 
		   lIsPosInOmegaProton && 
		   lIsNegInOmegaPion    )
		 	fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Omega Plus		
	if( lChargeXi > 0 && lAssoOmegaPlus) {
		lContainerPIDVars[0] = lmcPt              ;
		lContainerPIDVars[1] = lInvMassOmegaPlus  ;
		lContainerPIDVars[2] = lmcRapCasc         ;
		lContainerPIDVars[3] = lcentrality;
			
		// No PID
			fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorKaonForTPC  )
			fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
		if( lIsBachelorKaonForTPC && 
		    lIsNegProtonForTPC     )
			fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorKaonForTPC && 
		    lIsNegProtonForTPC    && 
		    lIsPosPionForTPC       )
			fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
		if( lIsBachelorKaon        )
			fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
		if( lIsBachelorKaon       && 
		    lIsNegInOmegaProton    )
			fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
		if( lIsBachelorKaon     && 
		    lIsNegInOmegaProton && 
		    lIsPosInOmegaPion    )
		 	fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
        
      }	
	
}// End of loop over reconstructed cascades
 
fHistnAssoXiMinus->Fill(nAssoXiMinus);
fHistnAssoXiPlus->Fill(nAssoXiPlus);
fHistnAssoOmegaMinus->Fill(nAssoOmegaMinus);
fHistnAssoOmegaPlus->Fill(nAssoOmegaPlus);  
//Printf("N asso Xi- = %n ", nAssoXiMinus);  
 
  // Post output data.
 PostData(1, fListHistCascade);
 PostData(2, fCFContCascadePIDAsXiMinus);
 PostData(3, fCFContCascadePIDAsXiPlus);
 PostData(4, fCFContCascadePIDAsOmegaMinus);
 PostData(5, fCFContCascadePIDAsOmegaPlus);
 PostData(6, fCFContAsCascadeCuts);

}      


//________________________________________________________________________
void AliAnalysisTaskCheckPerformanceCascadePbPb::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
	
  TList *cRetrievedList = 0x0;
  cRetrievedList = (TList*)GetOutputData(1);
  if(!cRetrievedList) {
	Printf("ERROR - AliAnalysisTaskCheckPerformanceCascadePbPb : ouput data container list not available\n");
	return;
  }	
	
  fHistMCTrackMultiplicity = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistMCTrackMultiplicity")  );
  if (!fHistMCTrackMultiplicity) {
    Printf("ERROR - AliAnalysisTaskCheckPerformanceCascadePbPb : fHistMCTrackMultiplicity not available");
    return;
  }
  
   
  TCanvas *canCheckPerformanceCascade = new TCanvas("CheckPerformanceCascadePbPb","Multiplicity",10,10,510,510);
  canCheckPerformanceCascade->cd(1)->SetLogy();

  fHistMCTrackMultiplicity->SetMarkerStyle(22);
  fHistMCTrackMultiplicity->DrawCopy("E");

}
