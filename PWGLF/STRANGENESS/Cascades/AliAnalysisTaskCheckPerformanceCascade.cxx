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
//		 AliAnalysisTaskCheckPerformanceCascade class
//            This task is for a performance study of cascade identification.
//            It works with MC info and ESD.
//              Use with AOD tree = under development
//            Origin   : A.Maire Mar2009, antonin.maire@ires.in2p3.fr
//            Modified : A.Maire Nov2010, antonin.maire@ires.in2p3.fr
//-----------------------------------------------------------------


#include <Riostream.h>

#include "TList.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
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
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
//   #include "AliV0vertexer.h"
//   #include "AliCascadeVertexer.h"
#include "AliESDEvent.h"
#include "AliESDcascade.h"

#include "AliAODEvent.h"

#include "AliAnalysisTaskCheckPerformanceCascade.h"

ClassImp(AliAnalysisTaskCheckPerformanceCascade)



     //_____Dummy constructor________________________________________________________________
AliAnalysisTaskCheckPerformanceCascade::AliAnalysisTaskCheckPerformanceCascade() 
: AliAnalysisTaskSE(), // <- take care to AliAnalysisTask( empty )
  fDebugCascade(0), fAnalysisType("ESD"), fTriggerMaskType("kMB"), fCollidingSystems(0), fESDpid(0), fESDtrackCuts(0), /*fPaveTextBookKeeping(0),*/
    fkRerunV0CascVertexers         (0),
    fkQualityCutZprimVtxPos        (kTRUE),
    fkRejectEventPileUp            (kTRUE),
    fkQualityCutNoTPConlyPrimVtx   (kTRUE),
    fkQualityCutTPCrefit           (kTRUE),
    fkQualityCut80TPCcls           (kTRUE),
    fkIsDataRecoWith1PadTPCCluster (kTRUE),
    fkExtraSelections              (0),
    
	// - Cascade part initialisation
fListHistCascade(0),
    fHistMCTrackMultiplicity(0), 
       // - Resolution of the multiplicity estimator
    f2dHistRecoPrimTrckMultVsMCMult(0), f2dHistRecoEstimateMultVsMCMult(0),

    fHistEtaGenProton(0),
    fHistEtaGenAntiProton(0),
    
   // Xi-
   fHistEtaGenCascXiMinus(0),
   f2dHistGenPtVsGenYGenXiMinus(0),
   
    fHistThetaGenCascXiMinus(0),
    f2dHistGenPtVsGenYFdblXiMinus(0),
    
    fHistThetaLambdaXiMinus(0), 
    fHistThetaBachXiMinus(0),
    
    fHistThetaMesDghterXiMinus(0), 
    fHistThetaBarDghterXiMinus(0),
    
    fHistPtBachXiMinus(0),
    fHistPtMesDghterXiMinus(0),
    fHistPtBarDghterXiMinus(0),
   
   
   // Xi+
   fHistEtaGenCascXiPlus(0),
   f2dHistGenPtVsGenYGenXiPlus(0),
   
    fHistThetaGenCascXiPlus(0), 
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
   f2dHistGenPtVsGenYGenOmegaMinus(0),
   
    fHistThetaGenCascOmegaMinus(0),
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
   f2dHistGenPtVsGenYGenOmegaPlus(0),
   
    fHistThetaGenCascOmegaPlus(0),
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

        // - Correlation in Pt between the cascade and its (anti)proton daughter
    f2dHistAsMCPtProtonVsPtXiMinus(0),
    f2dHistAsMCPtAntiProtonVsPtXiPlus(0),
    f2dHistAsMCPtProtonVsPtOmegaMinus(0),
    f2dHistAsMCPtAntiProtonVsPtOmegaPlus(0),

    fCFContCascadePIDAsXiMinus(0),
    fCFContCascadePIDAsXiPlus(0),
    fCFContCascadePIDAsOmegaMinus(0),
    fCFContCascadePIDAsOmegaPlus(0),

    fCFContAsCascadeCuts(0)

{
// Dummy constructor
        for(Int_t iAlephIdx   = 0; iAlephIdx   < 5; iAlephIdx++   ) { fAlephParameters [iAlephIdx]    = -1.; }
        for(Int_t iV0selIdx   = 0; iV0selIdx   < 7; iV0selIdx++   ) { fV0Sels          [iV0selIdx   ] = -1.; }
        for(Int_t iCascSelIdx = 0; iCascSelIdx < 8; iCascSelIdx++ ) { fCascSels        [iCascSelIdx ] = -1.; }
}
     
       
     
     
//_____Non-default Constructor________________________________________________________________
AliAnalysisTaskCheckPerformanceCascade::AliAnalysisTaskCheckPerformanceCascade(const char *name) 
  : AliAnalysisTaskSE(name),
    fDebugCascade(0), fAnalysisType("ESD"), fTriggerMaskType("kMB"), fCollidingSystems(0), fESDpid(0), fESDtrackCuts(0), /*fPaveTextBookKeeping(0),*/
    fkRerunV0CascVertexers         (0),
    fkQualityCutZprimVtxPos        (kTRUE),
    fkRejectEventPileUp            (kTRUE),
    fkQualityCutNoTPConlyPrimVtx   (kTRUE),
    fkQualityCutTPCrefit           (kTRUE),
    fkQualityCut80TPCcls           (kTRUE),
    fkIsDataRecoWith1PadTPCCluster (kTRUE),
    fkExtraSelections              (0),
      
    	// - Cascade part initialisation
fListHistCascade(0),
    fHistMCTrackMultiplicity(0), 
       // - Resolution of the multiplicity estimator
    f2dHistRecoPrimTrckMultVsMCMult(0), f2dHistRecoEstimateMultVsMCMult(0),

    fHistEtaGenProton(0),
    fHistEtaGenAntiProton(0),
    
// Xi-
   fHistEtaGenCascXiMinus(0),
   f2dHistGenPtVsGenYGenXiMinus(0),
   
    fHistThetaGenCascXiMinus(0),
    f2dHistGenPtVsGenYFdblXiMinus(0),
    
    fHistThetaLambdaXiMinus(0), 
    fHistThetaBachXiMinus(0),
    
    fHistThetaMesDghterXiMinus(0), 
    fHistThetaBarDghterXiMinus(0),
    
    fHistPtBachXiMinus(0),
    fHistPtMesDghterXiMinus(0),
    fHistPtBarDghterXiMinus(0),
   
   
   // Xi+
   fHistEtaGenCascXiPlus(0),
   f2dHistGenPtVsGenYGenXiPlus(0),
   
    fHistThetaGenCascXiPlus(0), 
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
   f2dHistGenPtVsGenYGenOmegaMinus(0),
   
    fHistThetaGenCascOmegaMinus(0),
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
   f2dHistGenPtVsGenYGenOmegaPlus(0),
   
    fHistThetaGenCascOmegaPlus(0),
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

        // - Correlation in Pt between the cascade and its (anti)proton daughter
    f2dHistAsMCPtProtonVsPtXiMinus(0),
    f2dHistAsMCPtAntiProtonVsPtXiPlus(0),
    f2dHistAsMCPtProtonVsPtOmegaMinus(0),
    f2dHistAsMCPtAntiProtonVsPtOmegaPlus(0),

    fCFContCascadePIDAsXiMinus(0),
    fCFContCascadePIDAsXiPlus(0),
    fCFContCascadePIDAsOmegaMinus(0),
    fCFContCascadePIDAsOmegaPlus(0),

    fCFContAsCascadeCuts(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // Output slot #1 writes into a TList container (cascade)
        
        for(Int_t iAlephIdx   = 0; iAlephIdx   < 5; iAlephIdx++   ) { fAlephParameters [iAlephIdx]    = -1.; }
        
        // New Loose : 1st step for the 7 TeV pp analysis
        
        fV0Sels[0] =  33.  ;  // max allowed chi2
        fV0Sels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
        fV0Sels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
        fV0Sels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
        fV0Sels[4] =   0.95;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
        fV0Sels[5] =   1.0 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
        fV0Sels[6] = 100.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
        
        fCascSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
        fCascSels[1] =   0.05 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
        fCascSels[2] =   0.010;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
        fCascSels[3] =   0.03 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
        fCascSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
        fCascSels[5] =   0.95 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
        fCascSels[6] =   0.4  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
        fCascSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
        
        
        // Hyper Loose "à la 900 GeV 2009 data", with lower cosine of pointing angle for Xi (0.95 down to 0.82) = 900 GeV paper
        /*
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
        */
        
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
        
  DefineOutput(1, TList::Class());
 
}


AliAnalysisTaskCheckPerformanceCascade::~AliAnalysisTaskCheckPerformanceCascade()
{
  //
  // Destructor
  //

  // For all TH1, 2, 3 HnSparse and CFContainer are in the fListCascade TList.
  // They will be deleted when fListCascade is deleted by the TSelector dtor
  // Because of TList::SetOwner()

  if (fListHistCascade)      { delete fListHistCascade;     fListHistCascade = 0x0;  }  
  if (fESDpid)               { delete fESDpid;              fESDpid = 0x0;} // fESDpid is not stored into the TList
  if (fESDtrackCuts)         { delete fESDtrackCuts;        fESDtrackCuts = 0x0; }
  /*if (fPaveTextBookKeeping)  { delete fPaveTextBookKeeping; fPaveTextBookKeeping = 0x0; } // fPaveTextBookKeeping is not stored into the TList*/
}


//________________________________________________________________________
void AliAnalysisTaskCheckPerformanceCascade::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

	
   // Option for AliLog
	AliLog::SetGlobalLogLevel(AliLog::kError); 
   	// to suppress the extensive info prompted by a run with MC			

   // Definition of the output datamembers	
   fListHistCascade = new TList();
   fListHistCascade->SetOwner(); // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
        
if(! fESDpid){

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
                 // fAlephParameters[0] = 4.23232575531564326e+00/50;//50*0.76176e-1;
                 // fAlephParameters[1] = 8.68482806165147636e+00;//10.632; 
                 // fAlephParameters[2] = 1.34000000000000005e-05;//0.13279e-4;
                 // fAlephParameters[3] = 2.30445734159456084e+00;//1.8631;
                 // fAlephParameters[4] = 2.25624744086878559e+00;//1.9479;
                
                // Param for LHC09d10 prod - A.Kalweit
                fAlephParameters[0] = 2.15898e+00/50.;
                fAlephParameters[1] = 1.75295e+01;
                fAlephParameters[2] = 3.40030e-09;
                fAlephParameters[3] = 1.96178e+00;
                fAlephParameters[4] = 3.91720e+00; 
         }
         Printf("CheckPerfCascade - Check Aleph Param in case of MC Data   (fAlephParameters[3] = %f) (To be compared with : 2.15 for 1-pad-cluster prod. / 1.96178 otherwise)\n",  fAlephParameters[3]);
        
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
        
        fListHistCascade->Add(fPaveTextBookKeeping);
}       
*/
                
  // - General
  
  if (!fHistMCTrackMultiplicity) {
     fHistMCTrackMultiplicity = new TH1F("fHistMCTrackMultiplicity", "MC Track Multiplicity;Number of MC tracks;Events", 100, 0, 500);
   //  fHistMCTrackMultiplicity = new TH1F("fHistMCTrackMultiplicity", "Multiplicity distribution;Number of tracks;Events", 200, 0, 40000); //HERE
  fListHistCascade->Add(fHistMCTrackMultiplicity);
  }
  
    // - Resolution of the multiplicity estimator
  if(! f2dHistRecoPrimTrckMultVsMCMult){
       f2dHistRecoPrimTrckMultVsMCMult = new TH2F("f2dHistRecoPrimTrckMultVsMCMult", "Resolution of the multiplicity estimator (prim. tracks in |#eta| < 0.8); Reco Multiplicity (prim. tracks); MC multiplicity (gen. part. in |#eta| < 0.8)", 120, 0., 120., 300, 0., 300.);
       fListHistCascade->Add(f2dHistRecoPrimTrckMultVsMCMult);
  }
  
  if(! f2dHistRecoEstimateMultVsMCMult){
       f2dHistRecoEstimateMultVsMCMult = new TH2F("f2dHistRecoEstimateMultVsMCMult", "Resolution of the multiplicity estimator (EstimateMult. in |#eta| < 1.0); Reco Multiplicity (tr(ITS-TPC)+ITSsa+tracklets); MC multiplicity (gen. part. in |#eta| < 1.0)", 160, 0., 160., 300, 0., 300.);
       fListHistCascade->Add(f2dHistRecoEstimateMultVsMCMult);
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
  
  if (!f2dHistGenPtVsGenYGenXiMinus) {
     f2dHistGenPtVsGenYGenXiMinus = new TH2F("f2dHistGenPtVsGenYGenXiMinus", "MC P_{t} Vs MC Y of Gen #Xi^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistGenPtVsGenYGenXiMinus);
  }
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascXiMinus) {
     fHistThetaGenCascXiMinus = new TH1F("fHistThetaGenCascXiMinus", "#theta of gen. #Xi^{-};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascXiMinus);
  }

  if (!f2dHistGenPtVsGenYFdblXiMinus) {
     f2dHistGenPtVsGenYFdblXiMinus = new TH2F("f2dHistGenPtVsGenYFdblXiMinus", "MC P_{t} Vs MC Y of findable Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
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
  
  
  
  //--------
  // II - Xi+ 
  // - Pseudo-Rapidity distribution
  if (!fHistEtaGenCascXiPlus) {
     fHistEtaGenCascXiPlus = new TH1F("fHistEtaGenCascXiPlus", "#eta of any gen. #bar{#Xi}^{+};#eta;Number of Casc", 200, -10, 10);
     fListHistCascade->Add(fHistEtaGenCascXiPlus);
  }
  
  if (!f2dHistGenPtVsGenYGenXiPlus) {
     f2dHistGenPtVsGenYGenXiPlus = new TH2F("f2dHistGenPtVsGenYGenXiPlus", "MC P_{t} Vs MC Y of Gen #bar{#Xi}^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistGenPtVsGenYGenXiPlus);
  }
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascXiPlus) {
     fHistThetaGenCascXiPlus = new TH1F("fHistThetaGenCascXiPlus", "#theta of gen. #bar{#Xi}^{+};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascXiPlus);
  }
 
  if (!f2dHistGenPtVsGenYFdblXiPlus) {
     f2dHistGenPtVsGenYFdblXiPlus = new TH2F("f2dHistGenPtVsGenYFdblXiPlus", "MC P_{t} Vs MC Y of findable Gen #bar{#Xi}^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
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
  
  if (!f2dHistGenPtVsGenYGenOmegaMinus) {
     f2dHistGenPtVsGenYGenOmegaMinus = new TH2F("f2dHistGenPtVsGenYGenOmegaMinus", "MC P_{t} Vs MC Y of Gen #Omega^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistGenPtVsGenYGenOmegaMinus);
  }
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascOmegaMinus) {
     fHistThetaGenCascOmegaMinus = new TH1F("fHistThetaGenCascOmegaMinus", "#theta of gen. #Omega^{-};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascOmegaMinus);
  }
 
  if (!f2dHistGenPtVsGenYFdblOmegaMinus) {
     f2dHistGenPtVsGenYFdblOmegaMinus = new TH2F("f2dHistGenPtVsGenYFdblOmegaMinus", "MC P_{t} Vs MC Y of findable Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
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
     fHistEtaGenCascOmegaPlus = new TH1F("fHistEtaGenCascOmegaPlus", "#eta of any gen. #bar{#Omega}^{+};#eta;Number of Casc", 200, -10, 10);
     fListHistCascade->Add(fHistEtaGenCascOmegaPlus);
  }
  
  if (!f2dHistGenPtVsGenYGenOmegaPlus) {
     f2dHistGenPtVsGenYGenOmegaPlus = new TH2F("f2dHistGenPtVsGenYGenOmegaPlus", "MC P_{t} Vs MC Y of Gen #bar{#Omega}^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistGenPtVsGenYGenOmegaPlus);
  }
  
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascOmegaPlus) {
     fHistThetaGenCascOmegaPlus = new TH1F("fHistThetaGenCascOmegaPlus", "#theta of gen. #bar{#Omega}^{+};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascOmegaPlus);
  }
 
  if (!f2dHistGenPtVsGenYFdblOmegaPlus) {
     f2dHistGenPtVsGenYFdblOmegaPlus = new TH2F("f2dHistGenPtVsGenYFdblOmegaPlus", "MC P_{t} Vs MC Y of findable Gen #bar{#Omega}^{+}; Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
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
	  fHistMassXiPlus = new TH1F("fHistMassXiPlus","#bar{#Xi}^{+} candidates;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",400,1.2,2.0);
	  fListHistCascade->Add(fHistMassXiPlus);
  }

  if (! fHistMassOmegaMinus) {
	  fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus","#Omega^{-} candidates;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500,1.5,2.5);
	  fListHistCascade->Add(fHistMassOmegaMinus);
  }
 
  if (! fHistMassOmegaPlus) {
	  fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus","#bar{#Omega}^{+} candidates;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 500,1.5,2.5);
	  fListHistCascade->Add(fHistMassOmegaPlus);
  }
  
  
  
  		// - Effective mass histos with combined PID
  
  if (! fHistMassWithCombPIDXiMinus) {
    fHistMassWithCombPIDXiMinus = new TH1F("fHistMassWithCombPIDXiMinus","#Xi^{-} candidates, with Bach. comb. PID;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 400,1.2,2.0);
    fListHistCascade->Add(fHistMassWithCombPIDXiMinus);
  }
  
  if (! fHistMassWithCombPIDXiPlus) {
    fHistMassWithCombPIDXiPlus = new TH1F("fHistMassWithCombPIDXiPlus","#bar{#Xi}^{+} candidates, with Bach. comb. PID;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",400,1.2,2.0);
    fListHistCascade->Add(fHistMassWithCombPIDXiPlus);
  }

  if (! fHistMassWithCombPIDOmegaMinus) {
	fHistMassWithCombPIDOmegaMinus = new TH1F("fHistMassWithCombPIDOmegaMinus","#Omega^{-} candidates, with Bach. comb. PID;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassWithCombPIDOmegaMinus);
  }
 
  if (! fHistMassWithCombPIDOmegaPlus) {
	fHistMassWithCombPIDOmegaPlus = new TH1F("fHistMassWithCombPIDOmegaPlus","#bar{#Omega}^{+} candidates, with Bach. comb. PID;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 500,1.5,2.5);
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
    fHistMassWithMcPIDXiPlus = new TH1F("fHistMassWithMcPIDXiPlus","#bar{#Xi}^{+} candidates, with Bach. MC PID;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",400,1.2,2.0);
    fListHistCascade->Add(fHistMassWithMcPIDXiPlus);
  }

  if (! fHistMassWithMcPIDOmegaMinus) {
	fHistMassWithMcPIDOmegaMinus = new TH1F("fHistMassWithMcPIDOmegaMinus","#Omega^{-} candidates, with Bach. MC PID;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassWithMcPIDOmegaMinus);
  }
 
  if (! fHistMassWithMcPIDOmegaPlus) {
	fHistMassWithMcPIDOmegaPlus = new TH1F("fHistMassWithMcPIDOmegaPlus","#bar{#Omega}^{+} candidates, with Bach. MC PID;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassWithMcPIDOmegaPlus);
  }
  
  
		// - Effective mass histos for cascades candidates ASSOCIATED with MC.
  
  if (! fHistAsMCMassXiMinus) {
	  fHistAsMCMassXiMinus = new TH1F("fHistAsMCMassXiMinus","#Xi^{-} candidates associated to MC;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 400,1.2,2.0);
	  fListHistCascade->Add(fHistAsMCMassXiMinus);
  }
  
  if (! fHistAsMCMassXiPlus) {
	  fHistAsMCMassXiPlus = new TH1F("fHistAsMCMassXiPlus","#bar{#Xi}^{+} candidates associated to MC;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",400,1.2,2.0);
	  fListHistCascade->Add(fHistAsMCMassXiPlus);
  }

  if (! fHistAsMCMassOmegaMinus) {
	  fHistAsMCMassOmegaMinus = new TH1F("fHistAsMCMassOmegaMinus","#Omega^{-} candidates associated to MC;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500,1.5,2.5);
	  fListHistCascade->Add(fHistAsMCMassOmegaMinus);
  }
 
  if (! fHistAsMCMassOmegaPlus) {
	  fHistAsMCMassOmegaPlus = new TH1F("fHistAsMCMassOmegaPlus","#bar{#Omega}^{+} candidates associated to MC;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 500,1.5,2.5);
	  fListHistCascade->Add(fHistAsMCMassOmegaPlus);
  }
  
		
		// -  Generated Pt Vs generated Y of the cascade candidates associated with MC 
		//     + having the proper maximum proba of combined PID for the bachelor
  
  if (!f2dHistAsMCandCombPIDGenPtVsGenYXiMinus) {
     f2dHistAsMCandCombPIDGenPtVsGenYXiMinus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYXiMinus", "MC P_{t} Vs MC Y of #Xi^{-} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYXiMinus);
  }
  
  if (!f2dHistAsMCandCombPIDGenPtVsGenYXiPlus) {
     f2dHistAsMCandCombPIDGenPtVsGenYXiPlus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYXiPlus", "MC P_{t} Vs MC Y of #bar{#Xi}^{+} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYXiPlus);
  }
  
  if (!f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus) {
     f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus", "MC P_{t} Vs MC Y of #Omega^{-} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus);
  }
  
  if (!f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus) {
     f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus", "MC P_{t} Vs MC Y of #bar{#Omega}^{+} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 220, -1.1, 1.1);
     fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus);
  }
  

  		// - Generated Pt Vs Generated Y, for the cascade candidates associated with MC
  
  if (!f2dHistAsMCGenPtVsGenYXiMinus) {
  	f2dHistAsMCGenPtVsGenYXiMinus = new TH2F("f2dHistAsMCGenPtVsGenYXiMinus", "MC P_{t} Vs MC Y of gen. #Xi^{-} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYXiMinus );
  }
  
  if (!f2dHistAsMCGenPtVsGenYXiPlus) {
	  f2dHistAsMCGenPtVsGenYXiPlus = new TH2F("f2dHistAsMCGenPtVsGenYXiPlus", "MC P_{t} Vs MC Y of gen. #bar{#Xi}^{+} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYXiPlus );
  }
  
  if (!f2dHistAsMCGenPtVsGenYOmegaMinus) {
	  f2dHistAsMCGenPtVsGenYOmegaMinus = new TH2F("f2dHistAsMCGenPtVsGenYOmegaMinus", "MC P_{t} Vs MC Y of gen. #Omega^{-} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYOmegaMinus );
  }
  
  if (!f2dHistAsMCGenPtVsGenYOmegaPlus) {
	  f2dHistAsMCGenPtVsGenYOmegaPlus = new TH2F("f2dHistAsMCGenPtVsGenYOmegaPlus", "MC P_{t} Vs MC Y of gen. #bar{#Omega}^{+} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 220, -1.1, 1.1);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYOmegaPlus );
  } 
  
  
  		// - Generated Eta of the the cascade candidates associated with MC
  if (!fHistAsMCGenEtaXiMinus) {
	  fHistAsMCGenEtaXiMinus = new TH1F("fHistAsMCGenEtaXiMinus", "#eta of gen. #Xi^{-} (associated);#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add( fHistAsMCGenEtaXiMinus );
  }
  
  if (!fHistAsMCGenEtaXiPlus) {
	  fHistAsMCGenEtaXiPlus = new TH1F("fHistAsMCGenEtaXiPlus", "#eta of gen. #bar{#Xi}^{+} (associated);#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add( fHistAsMCGenEtaXiPlus );
  }
  
  if (!fHistAsMCGenEtaOmegaMinus) {
	  fHistAsMCGenEtaOmegaMinus = new TH1F("fHistAsMCGenEtaOmegaMinus", "#eta of gen. #Omega^{-} (associated);#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add( fHistAsMCGenEtaOmegaMinus );
  }
  
  if (!fHistAsMCGenEtaOmegaPlus) {
	  fHistAsMCGenEtaOmegaPlus = new TH1F("fHistAsMCGenEtaOmegaPlus", "#eta of gen. #bar{#Omega}^{+} (associated);#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add( fHistAsMCGenEtaOmegaPlus );
  }
  
  
  
  		// - Resolution in Pt as function of generated Pt
  
  if(! f2dHistAsMCResPtXiMinus) {
	  f2dHistAsMCResPtXiMinus = new TH2F( "f2dHistAsMCResPtXiMinus", "Resolution in Pt reconstruction for #Xi^{-}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtXiMinus);
  }
  
  if(! f2dHistAsMCResPtXiPlus) {
	  f2dHistAsMCResPtXiPlus = new TH2F( "f2dHistAsMCResPtXiPlus", "Resolution in Pt reconstruction for #bar{#Xi}^{+}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtXiPlus);
  }
  
  if(! f2dHistAsMCResPtOmegaMinus) {
	  f2dHistAsMCResPtOmegaMinus = new TH2F( "f2dHistAsMCResPtOmegaMinus", "Resolution in Pt reconstruction for #Omega^{-}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtOmegaMinus);
  }
  
  if(! f2dHistAsMCResPtOmegaPlus) {
	  f2dHistAsMCResPtOmegaPlus = new TH2F( "f2dHistAsMCResPtOmegaPlus", "Resolution in Pt reconstruction for #bar{#Omega}^{+}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC}", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtOmegaPlus);
  }
  
  		// - Resolution in R(2D) as function of generated R
  
  if(! f2dHistAsMCResRXiMinus) {
	  f2dHistAsMCResRXiMinus = new TH2F( "f2dHistAsMCResRXiMinus", "Resolution in transv. position for #Xi^{-}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResRXiMinus);
  }
  
  if(! f2dHistAsMCResRXiPlus) {
	  f2dHistAsMCResRXiPlus = new TH2F( "f2dHistAsMCResRXiPlus", "Resolution in transv. position for #bar{#Xi}^{+}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResRXiPlus);
  }
  
  if(! f2dHistAsMCResROmegaMinus) {
	  f2dHistAsMCResROmegaMinus = new TH2F( "f2dHistAsMCResROmegaMinus", "Resolution in transv. position for #Omega^{-}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResROmegaMinus);
  }
  
  if(! f2dHistAsMCResROmegaPlus) {
	  f2dHistAsMCResROmegaPlus = new TH2F( "f2dHistAsMCResROmegaPlus", "Resolution in transv. position for #bar{#Omega}^{+}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC}", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResROmegaPlus);
  }
  
                // - Resolution in phi as function of generated Pt
    
  if(! f2dHistAsMCResPhiXiMinus) {
          f2dHistAsMCResPhiXiMinus = new TH2F( "f2dHistAsMCResPhiXiMinus", "Resolution in #phi for #Xi^{-}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);
          fListHistCascade->Add(f2dHistAsMCResPhiXiMinus);
  }
  
  if(! f2dHistAsMCResPhiXiPlus) {
          f2dHistAsMCResPhiXiPlus = new TH2F( "f2dHistAsMCResPhiXiPlus", "Resolution in #phi for #bar{#Xi}^{+}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);
          fListHistCascade->Add(f2dHistAsMCResPhiXiPlus);
  }
  
  if(! f2dHistAsMCResPhiOmegaMinus) {
          f2dHistAsMCResPhiOmegaMinus = new TH2F( "f2dHistAsMCResPhiOmegaMinus", "Resolution in #phi for #Omega^{-}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);  
          fListHistCascade->Add(f2dHistAsMCResPhiOmegaMinus);
  }
  
  if(! f2dHistAsMCResPhiOmegaPlus) {
          f2dHistAsMCResPhiOmegaPlus = new TH2F( "f2dHistAsMCResPhiOmegaPlus", "Resolution in #phi for #bar{#Omega}^{+}; Pt_{MC} (GeV/c); #phi(MC) - #phi(reco)   (deg)", 200, 0., 10., 60, -30., 30.);
          fListHistCascade->Add(f2dHistAsMCResPhiOmegaPlus);
  }
  
  
                // - Correlation in Pt between the cascade and its (anti)proton daughter
  if(! f2dHistAsMCPtProtonVsPtXiMinus) {
        f2dHistAsMCPtProtonVsPtXiMinus = new TH2F( "f2dHistAsMCPtProtonVsPtXiMinus", "Correlation Pt(p) Vs Pt(#Xi^{-}), associated to MC; Pt_{MC}(p) (GeV/c); Pt_{MC}(#Xi^{-}) (GeV/c)", 100, 0., 10., 100, 0., 10.);
        fListHistCascade->Add(f2dHistAsMCPtProtonVsPtXiMinus);
  }
  
  if(! f2dHistAsMCPtAntiProtonVsPtXiPlus) {
        f2dHistAsMCPtAntiProtonVsPtXiPlus = new TH2F( "f2dHistAsMCPtAntiProtonVsPtXiPlus", "Correlation Pt(#bar{p}) Vs Pt(#bar{#Xi}^{+}), associated to MC; Pt_{MC}(#bar{p}) (GeV/c); Pt_{MC}(#bar{#Xi}^{+}) (GeV/c)", 100, 0., 10., 100, 0., 10.);
        fListHistCascade->Add(f2dHistAsMCPtAntiProtonVsPtXiPlus);
  }
  
  if(! f2dHistAsMCPtProtonVsPtOmegaMinus) {
        f2dHistAsMCPtProtonVsPtOmegaMinus = new TH2F( "f2dHistAsMCPtProtonVsPtOmegaMinus", "Correlation Pt(p) Vs Pt(#Omega^{-}), associated to MC; Pt_{MC}(p) (GeV/c); Pt_{MC}(#Omega^{-}) (GeV/c)", 100, 0., 10., 100, 0., 10.);
        fListHistCascade->Add(f2dHistAsMCPtProtonVsPtOmegaMinus);
  }
  
  
  if(! f2dHistAsMCPtAntiProtonVsPtOmegaPlus) {
        f2dHistAsMCPtAntiProtonVsPtOmegaPlus = new TH2F( "f2dHistAsMCPtAntiProtonVsPtOmegaPlus", "Correlation Pt(#bar{p}) Vs Pt(#bar{#Omega}^{+}), associated to MC; Pt_{MC}(#bar{p}) (GeV/c); Pt_{MC}(#bar{#Omega}^{+}) (GeV/c)", 100, 0., 10., 100, 0., 10.);
        fListHistCascade->Add(f2dHistAsMCPtAntiProtonVsPtOmegaPlus);
  }
  
  
  
                // - PID container
if(! fCFContCascadePIDAsXiMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 200;
  lNbBinsPerVar[1] = 75;
  lNbBinsPerVar[2] = 44;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDAsXiMinus = new AliCFContainer("fCFContCascadePIDAsXiMinus","Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDAsXiMinus->SetBinLimits(0,   0.0  ,  10.0  );        // Pt(Cascade)
  fCFContCascadePIDAsXiMinus->SetBinLimits(1,   1.25 ,   1.40 );        // Xi Effective mass
  fCFContCascadePIDAsXiMinus->SetBinLimits(2,  -1.1  ,   1.1  );        // Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDAsXiMinus->SetBinLimits(3, 0.0, 20000.0  );    // nTrackPrimaryMultiplicity
  else
	fCFContCascadePIDAsXiMinus->SetBinLimits(3, 0.0, 250.0  );     // nTrackPrimaryMultiplicity
  
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
  fCFContCascadePIDAsXiMinus->SetVarTitle(3, "Primary Track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDAsXiMinus);
  
}

if(! fCFContCascadePIDAsXiPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 200;
  lNbBinsPerVar[1] = 75;
  lNbBinsPerVar[2] = 44;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDAsXiPlus = new AliCFContainer("fCFContCascadePIDAsXiPlus","Pt_{cascade} Vs M_{#bar{#Xi}^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDAsXiPlus->SetBinLimits(0,   0.0  ,  10.0  );         // Pt(Cascade)
  fCFContCascadePIDAsXiPlus->SetBinLimits(1,   1.25 ,   1.40 );         // Xi Effective mass
  fCFContCascadePIDAsXiPlus->SetBinLimits(2,  -1.1  ,   1.1  );         // Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDAsXiPlus->SetBinLimits(3, 0.0, 20000.0  );    // nTrackPrimaryMultiplicity
  else
	fCFContCascadePIDAsXiPlus->SetBinLimits(3, 0.0, 250.0  );     // nTrackPrimaryMultiplicity
  
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
  fCFContCascadePIDAsXiPlus->SetVarTitle(3, "Primary Track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDAsXiPlus);
  
}


if(! fCFContCascadePIDAsOmegaMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 200;
  lNbBinsPerVar[1] = 60;
  lNbBinsPerVar[2] = 44;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDAsOmegaMinus = new AliCFContainer("fCFContCascadePIDAsOmegaMinus","Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDAsOmegaMinus->SetBinLimits(0,   0.0  ,  10.0  );     // Pt(Cascade)
  fCFContCascadePIDAsOmegaMinus->SetBinLimits(1,   1.62 ,   1.74 );     // Omega Effective mass
  fCFContCascadePIDAsOmegaMinus->SetBinLimits(2,  -1.1  ,   1.1  );     // Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDAsOmegaMinus->SetBinLimits(3, 0.0, 20000.0  );    // nTrackPrimaryMultiplicity
  else
	fCFContCascadePIDAsOmegaMinus->SetBinLimits(3, 0.0, 250.0  );     // nTrackPrimaryMultiplicity
  
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
  fCFContCascadePIDAsOmegaMinus->SetVarTitle(3, "Primary Track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDAsOmegaMinus);
  
}

if(! fCFContCascadePIDAsOmegaPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4]= {0};
  lNbBinsPerVar[0] = 200;
  lNbBinsPerVar[1] = 60;
  lNbBinsPerVar[2] = 44;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDAsOmegaPlus = new AliCFContainer("fCFContCascadePIDAsOmegaPlus","Pt_{cascade} Vs M_{#bar{#Omega}^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDAsOmegaPlus->SetBinLimits(0,   0.0  ,  10.0  );      // Pt(Cascade)
  fCFContCascadePIDAsOmegaPlus->SetBinLimits(1,   1.62 ,   1.74 );      // Omega Effective mass
  fCFContCascadePIDAsOmegaPlus->SetBinLimits(2,  -1.1  ,   1.1  );      // Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDAsOmegaPlus->SetBinLimits(3, 0.0, 20000.0  );    // nTrackPrimaryMultiplicity
  else
	fCFContCascadePIDAsOmegaPlus->SetBinLimits(3, 0.0, 250.0  );     // nTrackPrimaryMultiplicity
  
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
  fCFContCascadePIDAsOmegaPlus->SetVarTitle(3, "Primary Track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDAsOmegaPlus);
  
}

// Part 3 : Towards the optimisation of topological selections -------
if(! fCFContAsCascadeCuts){
   
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
  
  lNbBinsPerVar[11] = 75; // 2-MeV/c2 bins
  lNbBinsPerVar[12] = 60; // 2-MeV/c2 bins
  
  lNbBinsPerVar[13] = 100;
  
  lNbBinsPerVar[14] = 44; // 0.05 in rapidity units
  lNbBinsPerVar[15] = 44; // 0.05 in rapidity units
  
  lNbBinsPerVar[16] = 20;
 
  lNbBinsPerVar[17] = 50;
  lNbBinsPerVar[18] = 100;
  lNbBinsPerVar[19] = 24;
   
   
  fCFContAsCascadeCuts = new AliCFContainer("fCFContAsCascadeCuts","Cut Container for Asso. Cascades", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  //0
  Double_t *lBinLim0  = new Double_t[ lNbBinsPerVar[0]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[0];i++)  lBinLim0[i]  = (Double_t)0.0   + (4.8  - 0.0 )/(lNbBinsPerVar[0]-1)  * (Double_t)i ;
        lBinLim0[ lNbBinsPerVar[0]  ] = 20.0;
  fCFContAsCascadeCuts -> SetBinLimits(0,  lBinLim0 );            // DcaXiDaughters : 0.0 to 5.0  
  delete [] lBinLim0;  
  //1
  Double_t *lBinLim1  = new Double_t[ lNbBinsPerVar[1]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[1];i++)   lBinLim1[i]  = (Double_t)0.0   + (0.24  - 0.0 )/(lNbBinsPerVar[1]-1)  * (Double_t)i ;
        lBinLim1[ lNbBinsPerVar[1]  ] = 100.0;
  fCFContAsCascadeCuts -> SetBinLimits(1,  lBinLim1 );            // DcaBachToPrimVertexXi : 0.0 to 0.25  
  delete [] lBinLim1;  
  //2 
  Double_t *lBinLim2  = new Double_t[ lNbBinsPerVar[2]+1 ];
        for(Int_t i=1; i< lNbBinsPerVar[2]+1;i++)   lBinLim2[i]  = (Double_t)0.81   + (1.0  - 0.81 )/(lNbBinsPerVar[2]-1)  * (Double_t) (i-1) ;   
        lBinLim2[0] = 0.0;
  fCFContAsCascadeCuts -> SetBinLimits(2,  lBinLim2 );            // XiCosineOfPointingAngle : 0.80 to 1.0        
  delete [] lBinLim2;  
  //3
  Double_t *lBinLim3  = new Double_t[ lNbBinsPerVar[3]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[3];i++)   lBinLim3[i]  = (Double_t)0.0   + (3.9  - 0.0 )/(lNbBinsPerVar[3]-1)  * (Double_t)i ;
        lBinLim3[ lNbBinsPerVar[3]  ] = 110.0;
  fCFContAsCascadeCuts -> SetBinLimits(3,  lBinLim3 );            // XiRadius : 0.0 to 4.0        
  delete [] lBinLim3;  
  //4
  fCFContAsCascadeCuts->SetBinLimits(4,    1.1  ,  1.13 );        // InvMassLambdaAsCascDghter
  //5
  Double_t *lBinLim5  = new Double_t[ lNbBinsPerVar[5]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[5];i++)   lBinLim5[i]  = (Double_t)0.0   + (4.8  - 0.0 )/(lNbBinsPerVar[5]-1)  * (Double_t)i ;
        lBinLim5[ lNbBinsPerVar[5]  ] = 20.0;
  fCFContAsCascadeCuts -> SetBinLimits(5,  lBinLim5 );            // DcaV0DaughtersXi : 0.0 to 5.0        
  delete [] lBinLim5;  
  
  
  //6
  Double_t *lBinLim6  = new Double_t[ lNbBinsPerVar[6]+1 ];
        for(Int_t i=1; i< lNbBinsPerVar[6]+1 ;i++)   lBinLim6[i]  = (Double_t)0.81   + (1.0  - 0.81 )/(lNbBinsPerVar[6]-1)  * (Double_t) (i-1) ;   
        lBinLim6[0] = 0.0;
  fCFContAsCascadeCuts -> SetBinLimits(6,  lBinLim6 );            // V0CosineOfPointingAngleXi : 0.80 to 1.0      
  delete [] lBinLim6;  
  //7
  Double_t *lBinLim7  = new Double_t[ lNbBinsPerVar[7]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[7];i++)   lBinLim7[i]  = (Double_t)0.0   + (7.8  - 0.0 )/(lNbBinsPerVar[7]-1)  * (Double_t)i ;
        lBinLim7[ lNbBinsPerVar[7]  ] = 100.0;
  fCFContAsCascadeCuts -> SetBinLimits(7,  lBinLim7 );            // V0RadiusXi : 0.0 to 8.0      
  delete [] lBinLim7;  
  //8
  Double_t *lBinLim8  = new Double_t[ lNbBinsPerVar[8]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[8];i++)   lBinLim8[i]  = (Double_t)0.0   + (0.39  - 0.0 )/(lNbBinsPerVar[8]-1)  * (Double_t)i ;
        lBinLim8[ lNbBinsPerVar[8]  ] = 100.0;
  fCFContAsCascadeCuts -> SetBinLimits(8,  lBinLim8 );            // DcaV0ToPrimVertexXi : 0.0 to 0.4     
  delete [] lBinLim8;  
  //9
  Double_t *lBinLim9  = new Double_t[ lNbBinsPerVar[9]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[9];i++)   lBinLim9[i]  = (Double_t)0.0   + (0.24  - 0.0 )/(lNbBinsPerVar[9]-1)  * (Double_t)i ;
        lBinLim9[ lNbBinsPerVar[9]  ] = 100.0;
  fCFContAsCascadeCuts -> SetBinLimits(9,  lBinLim9 );            // DcaPosToPrimVertexXi : 0.0 to 0.25   
  delete [] lBinLim9; 
  //10
  Double_t *lBinLim10  = new Double_t[ lNbBinsPerVar[10]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[10];i++)   lBinLim10[i]  = (Double_t)0.0   + (0.24  - 0.0 )/(lNbBinsPerVar[10]-1)  * (Double_t)i ;
        lBinLim10[ lNbBinsPerVar[10]  ] = 100.0;
  fCFContAsCascadeCuts -> SetBinLimits(10,  lBinLim10 );            // DcaPosToPrimVertexXi : 0.0 to 0.25 
  delete [] lBinLim10; 
  
  
  //11
  fCFContAsCascadeCuts->SetBinLimits(11,   1.25 ,  1.40 );	// InvMassXi
  fCFContAsCascadeCuts->SetBinLimits(12,   1.62  , 1.74 );	// InvMassOmega
  fCFContAsCascadeCuts->SetBinLimits(13,   0.0  , 10.0  );	// XiTransvMom 
  fCFContAsCascadeCuts->SetBinLimits(14,  -1.1  ,  1.1  );	// Y(Xi)
  fCFContAsCascadeCuts->SetBinLimits(15,  -1.1  ,  1.1  );	// Y(Omega)
  fCFContAsCascadeCuts->SetBinLimits(16, -10.0  , 10.0  );	// BestPrimaryVtxPosZ
  if(fCollidingSystems){
  	fCFContAsCascadeCuts->SetBinLimits(17,   0.0, 10000.0  );    // nTrackPrimaryMultiplicity
  	fCFContAsCascadeCuts->SetBinLimits(18,   0.0, 10000.0  );    // nITSandTPCtracksAndSPDtracklets
  }
  else{  
  	//fCFContAsCascadeCuts->SetBinLimits(17,   0.0, 250.0  );     // nTrackPrimaryMultiplicity
        Double_t *lBinLim17  = new Double_t[ lNbBinsPerVar[17]+1 ];
                lBinLim17[0] = 0;       lBinLim17[10] = 10;      lBinLim17[20] = 24;     lBinLim17[30] = 45;    lBinLim17[40] = 95;    lBinLim17[50] = 250;
                lBinLim17[1] = 1;       lBinLim17[11] = 11;      lBinLim17[21] = 25;     lBinLim17[31] = 50;    lBinLim17[41] = 100;
                lBinLim17[2] = 2;       lBinLim17[12] = 13;      lBinLim17[22] = 27;     lBinLim17[32] = 55;    lBinLim17[42] = 105;
                lBinLim17[3] = 3;       lBinLim17[13] = 14;      lBinLim17[23] = 30;     lBinLim17[33] = 60;    lBinLim17[43] = 110;
                lBinLim17[4] = 4;       lBinLim17[14] = 15;      lBinLim17[24] = 31;     lBinLim17[34] = 65;    lBinLim17[44] = 115;
                lBinLim17[5] = 5;       lBinLim17[15] = 16;      lBinLim17[25] = 32;     lBinLim17[35] = 70;    lBinLim17[45] = 120;
                lBinLim17[6] = 6;       lBinLim17[16] = 20;      lBinLim17[26] = 33;     lBinLim17[36] = 75;    lBinLim17[46] = 125;
                lBinLim17[7] = 7;       lBinLim17[17] = 21;      lBinLim17[27] = 34;     lBinLim17[37] = 80;    lBinLim17[47] = 130;
                lBinLim17[8] = 8;       lBinLim17[18] = 22;      lBinLim17[28] = 35;     lBinLim17[38] = 85;    lBinLim17[48] = 135;
                lBinLim17[9] = 9;       lBinLim17[19] = 23;      lBinLim17[29] = 40;     lBinLim17[39] = 90;    lBinLim17[49] = 140;
                
                fCFContAsCascadeCuts -> SetBinLimits(17,  lBinLim17 );            // nTrackPrimaryMultiplicity : 0 to 250	
        delete [] lBinLim17;     
          

  	fCFContAsCascadeCuts->SetBinLimits(18,   0.0, 200.0  );     // nITSandTPCtracksAndSPDtracklets
  }
  fCFContAsCascadeCuts->SetBinLimits(19,  68.0  ,164.0  );	// BachTPCClusters
  
  
  // Regular binning definition (valid for v4-18-10-AN on)
  /*
  //setting the bin limits
  fCFContAsCascadeCuts->SetBinLimits(0,    0.0  ,  2.5  );	// DcaXiDaughters
  fCFContAsCascadeCuts->SetBinLimits(1,    0.0  ,  0.25 );	// DcaBachToPrimVertexXi
  fCFContAsCascadeCuts->SetBinLimits(2,    0.99 ,  1.0  );	// XiCosineOfPointingAngle
  fCFContAsCascadeCuts->SetBinLimits(3,    0.0  ,  4.0  );	// XiRadius
  fCFContAsCascadeCuts->SetBinLimits(4,    1.1  ,  1.15 );	// InvMassLambdaAsCascDghter
  fCFContAsCascadeCuts->SetBinLimits(5,    0.0  ,  1.0  );	// DcaV0DaughtersXi
  fCFContAsCascadeCuts->SetBinLimits(6,    0.98 ,  1.0  );	// V0CosineOfPointingAngleXi
  fCFContAsCascadeCuts->SetBinLimits(7,    0.0  , 20.0  );	// V0RadiusXi
  fCFContAsCascadeCuts->SetBinLimits(8,    0.0  ,  1.0  );	// DcaV0ToPrimVertexXi
  fCFContAsCascadeCuts->SetBinLimits(9,    0.0  ,  0.25 );	// DcaPosToPrimVertexXi
  fCFContAsCascadeCuts->SetBinLimits(10,   0.0  ,  0.25 );	// DcaNegToPrimVertexXi
  fCFContAsCascadeCuts->SetBinLimits(11,   1.25 ,  1.40 );	// InvMassXi
  fCFContAsCascadeCuts->SetBinLimits(12,   1.62 ,  1.74 );	// InvMassOmega
  fCFContAsCascadeCuts->SetBinLimits(13,   0.0  , 10.0  );	// pt_MC(Xi)
  fCFContAsCascadeCuts->SetBinLimits(14,  -1.1  ,  1.1  );	// Y_MC(Xi)
  fCFContAsCascadeCuts->SetBinLimits(15,  -1.1  ,  1.1  );	// Y_MC(Omega)
  fCFContAsCascadeCuts->SetBinLimits(16, -10.0  , 10.0  );	// BestPrimaryVtxPosZ
  if(fCollidingSystems){
  	fCFContAsCascadeCuts->SetBinLimits(17,   0.0, 10000.0  );    // nTrackPrimaryMultiplicity
  	fCFContAsCascadeCuts->SetBinLimits(18,   0.0, 10000.0  );    // nITSandTPCtracksAndSPDtracklets
  }
  else{  
  	fCFContAsCascadeCuts->SetBinLimits(17,   0.0, 250.0  );     // nTrackPrimaryMultiplicity
  	fCFContAsCascadeCuts->SetBinLimits(18,   0.0, 200.0  );     // nITSandTPCtracksAndSPDtracklets
  }
  fCFContAsCascadeCuts->SetBinLimits(19,  25.0  ,165.0  );	// BachTPCClusters
  */
  
  // Setting the number of steps : one for each cascade species (Xi-, Xi+ and Omega-, Omega+)
  fCFContAsCascadeCuts->SetStepTitle(0, "#Xi^{-} candidates associated to MC");
  fCFContAsCascadeCuts->SetStepTitle(1, "#bar{#Xi}^{+} candidates associated to MC");
  fCFContAsCascadeCuts->SetStepTitle(2, "#Omega^{-} candidates associated to MC");
  fCFContAsCascadeCuts->SetStepTitle(3, "#bar{#Omega}^{+} candidates associated to MC");
  
  // Setting the variable title, per axis
  // fCFContAsCascadeCuts->SetVarTitle(40,  "Chi2Xi");
  fCFContAsCascadeCuts->SetVarTitle(0,  "Dca(XiDaughters) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(1,  "Dca(Bach/PrimVertex) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(2,  "cos(Xi pointing angle)");
  fCFContAsCascadeCuts->SetVarTitle(3,  "R_{2d}(Xi decay) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(4,  "M_{#Lambda}(As Casc Dghter) (GeV/c^{2})");
   // fCFContAsCascadeCuts->SetVarTitle(40,  "V0Chi2Xi");
  fCFContAsCascadeCuts->SetVarTitle(5,  "Dca(V0 Daughters Xi) (cm)");
  
  fCFContAsCascadeCuts->SetVarTitle(6,  "cos(V0 pointing Angle) in Casc");
  fCFContAsCascadeCuts->SetVarTitle(7,  "R_{2d}(V0 decay) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(8,  "Dca(V0/PrimVertex) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(9,  "Dca(Pos/PrimVertex) (cm)");
  fCFContAsCascadeCuts->SetVarTitle(10, "Dca(Neg/PrimVertex) (cm)");
  
  fCFContAsCascadeCuts->SetVarTitle(11, "Inv. Mass(Xi) (GeV/c^{2})");
  fCFContAsCascadeCuts->SetVarTitle(12, "Inv. Mass(Omega) (GeV/c^{2})");
  
  fCFContAsCascadeCuts->SetVarTitle(13, "Pt_{MC}(Casc.) (GeV/c)");
  //fCFContAsCascadeCuts->SetVarTitle(40, "V0toXiCosineOfPointingAngle");
  
  fCFContAsCascadeCuts->SetVarTitle(14, "Y_{MC}(Xi)");
  fCFContAsCascadeCuts->SetVarTitle(15, "Y_{MC}(Omega)");
  
  fCFContAsCascadeCuts->SetVarTitle(16, "Z-position(BestPrimVtx) (cm)");
  
  fCFContAsCascadeCuts->SetVarTitle(17, "Primary Track Multiplicity");
  fCFContAsCascadeCuts->SetVarTitle(18, "(ITS+TPC tracks + SPD tracklets) Multiplicity");
  fCFContAsCascadeCuts->SetVarTitle(19, "Bach.TPC Clusters");
  
  fListHistCascade->Add(fCFContAsCascadeCuts);
}


PostData(1, fListHistCascade); 
}// end CreateOutputObjects






//________________________________________________________________________
void AliAnalysisTaskCheckPerformanceCascade::UserExec(Option_t *) 
{
	
  // Main loop
  // Called for each event
	
	AliESDEvent *lESDevent = 0x0;
	AliAODEvent *lAODevent = 0x0;
	AliMCEvent  *lMCevent  = 0x0; 
	AliStack    *lMCstack  = 0x0; 
	Int_t ncascades = -1;
	
	
  // Connect to the InputEvent	
  // After these lines, we should have an ESD/AOD event + the number of cascades in it.
		
	if(fAnalysisType == "ESD"){
		lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
		if (!lESDevent) {
			Printf("ERROR: lESDevent not available \n");
			cout << "Name of the file with pb :" <<  CurrentFileName() << endl;  // or AliAnalysisTaskSE::CurrentFileName()
			return;
		}
	}
  
	else if(fAnalysisType == "AOD"){  
		lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() ); 
		if (!lAODevent) {
			Printf("ERROR: lAODevent not available \n");
			cout << "Name of the file with pb :" <<  CurrentFileName() << endl;
			return;
		}
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
	
        
   // Temporary way : AOD awareness of the code to be developed  FIXME
   if(fAnalysisType == "AOD") return;
   
   
   
   
  //-------------------------------------------------
  // 0 - Trigger managment + global event selection
  // NOTE : Check the availability of the proper trigger 

   // Note : Presuppose the presence of AliPhysicsSelectionTask
   
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
        }


  //------------------------------------------------
  // 2 - Preparing the general info about of the event = prim. Vtx + magnetic field (ESD)
  

//      if(fAnalysisType == "ESD" ){

        // Magnetic field
                const Double_t lMagneticField = lESDevent->GetMagneticField( );

        // Prim vertex
                const AliESDVertex *lPrimaryTrackingVtx = lESDevent->GetPrimaryVertexTracks();  // get the vtx stored in ESD found with tracks
                const AliESDVertex *lPrimarySPDVtx      = lESDevent->GetPrimaryVertexSPD();     // get the vtx stored in ESD found with SPD tracklets

                const AliESDVertex *lPrimaryBestVtx     = lESDevent->GetPrimaryVertex();
                        // get the best primary vertex available for the event
                        // As done in AliCascadeVertexer, we keep the one which is the best one available.
                        // between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
                        Double_t lBestPrimaryVtxPos[3]   = {-100.0, -100.0, -100.0};
                lPrimaryBestVtx->GetXYZ( lBestPrimaryVtxPos );

                // FIXME : quality cut on the z-position of the prim vertex.
                if(fkQualityCutZprimVtxPos) {
                        if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0 ) { 
                                AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !"); 
                                PostData(1, fListHistCascade); 
                                return;
                        }
                }
                // FIXME : quality selection regarding pile-up rejection 
                if(fkRejectEventPileUp) {
                        if(lESDevent->IsPileupFromSPDInMultBins() ){// minContributors=3, minZdist=0.8, nSigmaZdist=3., nSigmaDiamXY=2., nSigmaDiamZ=5.  -> see http://alisoft.cern.ch/viewvc/trunk/STEER/AliESDEvent.h?root=AliRoot&r1=41914&r2=42199&pathrev=42199
                                AliWarning("Pb / Event tagged as pile-up by SPD... return !"); 
                                PostData(1, fListHistCascade); 
                                return; 
                        }
                }
                // FIXME : remove TPC-only primary vertex : retain only events with tracking + SPD vertex
                if(fkQualityCutNoTPConlyPrimVtx) {
                        if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingVtx->GetStatus() ){
                                AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
                                PostData(1, fListHistCascade); 
                                return;
                        }
                }
//      }// if ESD
        
        
  //	cout << "Name of the accessed file :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;

  //	cout << "Tree characteristics ..." << endl;
  //	fInputHandler->GetTree()->Print("toponly");
  //	fInputHandler->GetTree()->GetBranch("PrimaryVertex")->Print();
  //	fInputHandler->GetTree()->GetBranch("SPDVertex")->Print();



  // ---------------------------------------------------------------
  // - Initialisation of the part dedicated to cascade vertices

  if(fAnalysisType == "ESD")            ncascades = lESDevent->GetNumberOfCascades();
  else if(fAnalysisType == "AOD")       ncascades = lAODevent->GetNumberOfCascades();
	
  
  Int_t   nNumberOfMCPrimaries       = -1;
  Int_t   nMCPrimariesInEtaBelow0p8  =  0;
  Int_t   nMCPrimariesInEtaBelow1p0  =  0;
  
  Int_t   nTrackPrimaryMultiplicity        = -1;
  Int_t   nSPDTracklets                    =  0; // AliESDEvent::EstimateMultiplicity will re-initialise the value to 0
  Int_t   nITSandTPCtracksAndSPDtracklets  =  0; // AliESDEvent::EstimateMultiplicity will re-initialise the value to 0
  Int_t   nTracksITSSApure                 =  0; // AliESDEvent::EstimateMultiplicity will re-initialise the value to 0
  
        nNumberOfMCPrimaries       = lMCstack->GetNprimary();
        if(nNumberOfMCPrimaries < 1) return;
        
        nTrackPrimaryMultiplicity  = fESDtrackCuts->CountAcceptedTracks(lESDevent);
        //EstimateMultiplicity(Int_t &tracklets, Int_t &trITSTPC, Int_t &trITSSApure, Double_t eta, Bool_t useDCAFlag,Bool_t useV0Flag)
        lESDevent->EstimateMultiplicity( nSPDTracklets, nITSandTPCtracksAndSPDtracklets, nTracksITSSApure, 1.0, kTRUE, kTRUE);


        fHistMCTrackMultiplicity->Fill( nNumberOfMCPrimaries );
    
//_____________________________________________________________________________	
// Part 1 - Loop over the MC primaries	
        
    for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < nNumberOfMCPrimaries; iCurrentLabelStack++) 
    {// This is the begining of the loop on primaries
      
        TParticle* lCurrentParticle = 0x0; 
                   lCurrentParticle = lMCstack->Particle( iCurrentLabelStack );
        if(!lCurrentParticle){
                Printf("MC Primary loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
                continue;
        }
        
        Double_t lEtaCurrentParticle = TMath::Abs( lCurrentParticle->Eta() );
        if( lEtaCurrentParticle < 1.0 ){
                nMCPrimariesInEtaBelow1p0++;
                if( lEtaCurrentParticle < 0.8 ) nMCPrimariesInEtaBelow0p8++;
        }
    }
    
    f2dHistRecoPrimTrckMultVsMCMult->Fill( nTrackPrimaryMultiplicity,       nMCPrimariesInEtaBelow0p8 );
    f2dHistRecoEstimateMultVsMCMult->Fill( nITSandTPCtracksAndSPDtracklets, nMCPrimariesInEtaBelow1p0 );


   // For proton
   /*
   for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < nNumberOfMCPrimaries; iCurrentLabelStack++) 
    	{// This is the begining of the loop on primaries, for protons
          
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
	TH2F *l2dHistGenPtVsGenYGen   = 0;
		
	TH1F *lHistThetaGenCasc       = 0;
	TH2F *l2dHistGenPtVsGenYFdbl  = 0;
	TH1F *lHistThetaLambda        = 0;
	TH1F *lHistThetaBach          = 0;
	TH1F *lHistThetaBarDghter     = 0;
	TH1F *lHistThetaMesDghter     = 0;
	TH1F *lHistPtBach             = 0;
	TH1F *lHistPtBarDghter        = 0;
	TH1F *lHistPtMesDghter        = 0;


for(Int_t iCascType = 1; iCascType < 5; iCascType++)
{
       
switch (iCascType)
  {
    case 1: // Xi-
         lPdgCodeCasc       =   3312;  //Xi-
         lPdgCodeBach       =   -211;  //Pi-
         lPdgCodeLambda     =   3122;  //Lambda0
         lPdgCodeDghtMesV0  =   -211;  //Pi-
         lPdgCodeDghtBarV0  =   2212;  //Proton 
	 	
	 	// any Xi-
	 lHistEtaGenCasc        = fHistEtaGenCascXiMinus;
	 l2dHistGenPtVsGenYGen  = f2dHistGenPtVsGenYGenXiMinus;
	
	 	// cascades generated within acceptance (cut in pt + theta)
	 lHistThetaGenCasc      = fHistThetaGenCascXiMinus;
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
	 l2dHistGenPtVsGenYGen  = f2dHistGenPtVsGenYGenXiPlus;
	
	 	// cascades generated within acceptance (cut in pt + theta)	 
	 lHistThetaGenCasc      = fHistThetaGenCascXiPlus;
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
	 l2dHistGenPtVsGenYGen  = f2dHistGenPtVsGenYGenOmegaMinus;	
	
	 	// cascades generated within acceptance (cut in pt + theta)
	 lHistThetaGenCasc      = fHistThetaGenCascOmegaMinus;
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
	 l2dHistGenPtVsGenYGen  = f2dHistGenPtVsGenYGenOmegaPlus;		
	 	
	 	// cascades generated within acceptance (cut in pt + theta)
	 lHistThetaGenCasc      = fHistThetaGenCascOmegaPlus;
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


   for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < nNumberOfMCPrimaries; iCurrentLabelStack++) 
    {// This is the begining of the loop on primaries
      
        TParticle* lCurrentParticle = 0x0; 
    	           lCurrentParticle = lMCstack->Particle( iCurrentLabelStack );
	if(!lCurrentParticle){
		Printf("Cascade loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
		continue;
		
	}
	 
    	if( lCurrentParticle->GetPdgCode() == lPdgCodeCasc ){  // Here !
   		//cout << "Xi- within loop " << iCurrentLabelStack << "/ " << nNumberOfMCPrimaries << endl;
		
		// -  Xi level ... _____________________________________________________________
		TParticle* xiMC = 0x0;
			   xiMC = lCurrentParticle;
		if(!xiMC){
			Printf("MC TParticle pointer to Cascade = 0x0 ! Skip ...");
			continue;
		
		}
		
		// Fill the first histos : = any generated Xi, not necessarily within the acceptance
		Double_t lRapXiMC = 0.5*TMath::Log((xiMC->Energy() + xiMC->Pz()) / (xiMC->Energy() - xiMC->Pz() +1.e-13));
		
		lHistEtaGenCasc 	->Fill( xiMC->Eta() );	 
		l2dHistGenPtVsGenYGen 	->Fill( xiMC->Pt(), lRapXiMC  );    	
			
		
		
		// Check the emission of particle stays within the acceptance of the detector (cut in theta)
		if( xiMC->Theta() < TMath::Pi()/4.0  ||    xiMC->Theta() > 3.0*TMath::Pi()/4.0 ) continue;	
		if( xiMC->GetNDaughters() != 2) continue;
		if( xiMC->GetDaughter(0) < 0 )  continue;
		if( xiMC->GetDaughter(1) < 0 )  continue;
		
			TParticle* lDght0ofXi = lMCstack->Particle(  xiMC->GetDaughter(0) );
			TParticle* lDght1ofXi = lMCstack->Particle(  xiMC->GetDaughter(1) );
			
		TParticle* lLambda = 0;
		TParticle* lBach   = 0;
			
		// Xi - Case 1
			if(	lDght0ofXi->GetPdgCode() == lPdgCodeLambda   &&  // Here !
				lDght1ofXi->GetPdgCode() == lPdgCodeBach ){      // Here !
				
				lLambda = lDght0ofXi;
				lBach   = lDght1ofXi;
			}// end if dghter 0 = Lambda and    dghter 1 = Pi-  
			
		// Xi - Case 2
			else if( lDght0ofXi->GetPdgCode() == lPdgCodeBach  &&      // Here !
				 lDght1ofXi->GetPdgCode() == lPdgCodeLambda ){     // Here !
			
				lBach   = lDght0ofXi;
				lLambda = lDght1ofXi;
			}//  end if dghter 0 = Pi-  and   dghter 1 = Lambda
			
		// V0 otherwise - Case 3	
			else continue;
		
			// Check the emission of particle stays within the acceptance of the detector (cut in pt + theta)
			if( lLambda->Theta() < TMath::Pi()/4.0  ||    lLambda->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
			if( lBach->Theta() < TMath::Pi()/4.0    ||    lBach->Theta() > 3.0*TMath::Pi()/4.0 )   continue;
		
			if( lBach->Pt() < 0.150 ) continue; //FIXME : maybe tuned for Xi but not for K- from Omega ...
			
		
		
		// -  V0 level ... _____________________________________________________________
		TParticle* lDghtBarV0 = 0;
		TParticle* lDghtMesV0 = 0;
		
		if( lLambda->GetNDaughters() != 2 )  continue;
		if( lLambda->GetDaughter(0) < 0 )    continue;
		if( lLambda->GetDaughter(1) < 0 )    continue;
		
		
		TParticle* lDght0ofLambda = lMCstack->Particle(  lLambda->GetDaughter(0) );
		TParticle* lDght1ofLambda = lMCstack->Particle(  lLambda->GetDaughter(1) );
			
		// V0 - Case 1
			if( 	lDght0ofLambda->GetPdgCode() == lPdgCodeDghtBarV0 &&    // Here !
				lDght1ofLambda->GetPdgCode() == lPdgCodeDghtMesV0 ){    // Here !
			
				lDghtBarV0 = lDght0ofLambda;
				lDghtMesV0 = lDght1ofLambda;
			}// end if dghter 0 = Proton  and   dghter 1 = Pi-  
			
		// V0 - Case 2
			else if( lDght0ofLambda->GetPdgCode() == lPdgCodeDghtMesV0  &&     // Here !
				 lDght1ofLambda->GetPdgCode() == lPdgCodeDghtBarV0 ){      // Here !
			
				lDghtMesV0 = lDght0ofLambda;
				lDghtBarV0 = lDght1ofLambda;
			}//  end if dghter 0 = Pi-  and   dghter 1 = proton
			
		// V0 otherwise - Case 3
			else continue;
	
			
			// Check the emission of particle stays within the acceptance of the detector
			if( lDghtBarV0->Theta() < TMath::Pi()/4.0  ||  lDghtBarV0->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
			if( lDghtMesV0->Theta() < TMath::Pi()/4.0  ||  lDghtMesV0->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
			
			if( lDghtBarV0->Pt() < 0.250 ) continue;
			if( lDghtMesV0->Pt() < 0.150 ) continue;
			
			
			
		// - Just to know which file is currently open : locate the file containing Xi 
		//cout << "Name of the file containing generated Xi :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() 
		//						     <<  endl;	
			
		Double_t lRadToDeg = 180.0/TMath::Pi();	
			
		// - Filling histos ... _________________________________________________________________	
			lHistThetaGenCasc	->Fill( lRadToDeg * xiMC->Theta()  );
			l2dHistGenPtVsGenYFdbl	->Fill( xiMC->Pt(), lRapXiMC );
			
			// - Fill theta histos for Lambda and Bach
			lHistThetaLambda	->Fill( lRadToDeg * lLambda->Theta() );
			lHistThetaBach  	->Fill( lRadToDeg *   lBach->Theta() );
			
			// - Fill theta histos for V0 daughters
			lHistThetaBarDghter	->Fill( lRadToDeg * lDghtBarV0->Theta() );
			lHistThetaMesDghter	->Fill( lRadToDeg * lDghtMesV0->Theta() );
			
			// - Fill pt histos.
			lHistPtBach		->Fill(      lBach->Pt() );
			lHistPtBarDghter	->Fill( lDghtBarV0->Pt() );
			lHistPtMesDghter	->Fill( lDghtMesV0->Pt() );
						
		}// end if current particle = Xi-
	     
     }// This is the end of the loop on primaries
     
// - Re-initialisation of the local TH1F pointers
lHistEtaGenCasc         = 0x0;
l2dHistGenPtVsGenYGen   = 0x0;

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
 	
 
 
//__________________________________________________________________________	
// Part 3 - Loop over the reconstructed candidates
  

for (Int_t iXi = 0; iXi < ncascades; iXi++) 
{// This is the begining of the Cascade loop
		
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
	
        Int_t lPosTPCClusters   = pTrackXi->GetTPCNcls();
        Int_t lNegTPCClusters   = nTrackXi->GetTPCNcls();
        Int_t lBachTPCClusters  = bachTrackXi->GetTPCNcls(); 

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
	
	// - Step II.2 : Info over reconstructed cascades
	//-------------	
	
	Double_t lInvMassXiMinus    = 0.;
	Double_t lInvMassXiPlus     = 0.;
	Double_t lInvMassOmegaMinus = 0.;
	Double_t lInvMassOmegaPlus  = 0.;
	
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
		
	Double_t lChargeXi = xiESD->Charge();
	
	if( lChargeXi < 0 )	fHistMassXiMinus	->Fill( lInvMassXiMinus );
	if( lChargeXi > 0 )	fHistMassXiPlus		->Fill( lInvMassXiPlus );
	if( lChargeXi < 0 )	fHistMassOmegaMinus	->Fill( lInvMassOmegaMinus );
	if( lChargeXi > 0 )	fHistMassOmegaPlus	->Fill( lInvMassOmegaPlus );
	
	
	// - Step II.3 : PID info
	//-------------
	
	
	// 3.1 - PID Information

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

        // 3.1.A - Combined PID
	// Reasonable guess for the priors for the cascade track sample (e-, mu, pi, K, p)
	Double_t lPriorsGuessXi[10]    = {0, 0, 2, 0, 1,  0,0,0,0,0};
	Double_t lPriorsGuessOmega[10] = {0, 0, 1, 1, 1,  0,0,0,0,0};
	
	// Combined VO-positive-daughter PID
	AliPID pPidXi;         pPidXi.SetPriors(    lPriorsGuessXi    , kTRUE); // kTRUE = for charged particle PID
	AliPID pPidOmega;      pPidOmega.SetPriors( lPriorsGuessOmega , kTRUE); // kTRUE = for charged particle PID
		
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
	AliPID nPidXi;         nPidXi.SetPriors(    lPriorsGuessXi    , kTRUE); // kTRUE = for charged particle PID
	AliPID nPidOmega;      nPidOmega.SetPriors( lPriorsGuessOmega , kTRUE); // kTRUE = for charged particle PID
		
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
	AliPID bachPidXi;      bachPidXi.SetPriors(    lPriorsGuessXi    , kTRUE); // kTRUE = for charged particle PID
	AliPID bachPidOmega;   bachPidOmega.SetPriors( lPriorsGuessOmega , kTRUE); // kTRUE = for charged particle PID
	
        Double_t ppionBach = 0.0, pkaonBach = 0.0;
        
	if( bachTrackXi->IsOn(AliESDtrack::kESDpid) ){  // Combined PID exists
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
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 4) lIsBachelorKaonForTPC = kTRUE;
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 4) lIsBachelorPionForTPC = kTRUE;
        
        // Negative V0 daughter
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 4) lIsNegPionForTPC   = kTRUE;
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 4) lIsNegProtonForTPC = kTRUE;
        
        // Positive V0 daughter
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 4) lIsPosPionForTPC   = kTRUE;
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 4) lIsPosProtonForTPC = kTRUE;
        
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
		
        // Combined PID TH1s
	if( lChargeXi < 0 && lIsBachelorPion )    fHistMassWithCombPIDXiMinus     ->Fill( lInvMassXiMinus    );
	if( lChargeXi > 0 && lIsBachelorPion )    fHistMassWithCombPIDXiPlus      ->Fill( lInvMassXiPlus     );
	if( lChargeXi < 0 && lIsBachelorKaon )    fHistMassWithCombPIDOmegaMinus  ->Fill( lInvMassOmegaMinus );
	if( lChargeXi > 0 && lIsBachelorKaon )    fHistMassWithCombPIDOmegaPlus   ->Fill( lInvMassOmegaPlus  );
         
	
	// 3.2 - PID proba Vs Pt(Bach)
	Int_t      lblBachForPID  = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
	TParticle* mcBachForPID   = lMCstack->Particle( lblBachForPID );
	Double_t   lmcPtBach      = mcBachForPID->Pt();
	
	if(lIsBachelorPion)   f2dHistPIDprobaPionVsMCPtBach->Fill( lmcPtBach, ppionBach );
	if(lIsBachelorKaon)   f2dHistPIDprobaKaonVsMCPtBach->Fill( lmcPtBach, pkaonBach );
	
			
	// 3.3 - MC perfect PID
	Bool_t   lIsBachelorMCPiMinus  = kFALSE;
	Bool_t   lIsBachelorMCPiPlus   = kFALSE;
	Bool_t   lIsBachelorMCKMinus   = kFALSE;
	Bool_t   lIsBachelorMCKPlus    = kFALSE;	
	
	if( mcBachForPID->GetPdgCode() == -211) lIsBachelorMCPiMinus = kTRUE;
	if( mcBachForPID->GetPdgCode() ==  211) lIsBachelorMCPiPlus  = kTRUE;
	if( mcBachForPID->GetPdgCode() == -321) lIsBachelorMCKMinus  = kTRUE;
	if( mcBachForPID->GetPdgCode() ==  321) lIsBachelorMCKPlus   = kTRUE;
	
	if( lChargeXi < 0 && lIsBachelorMCPiMinus )    fHistMassWithMcPIDXiMinus     ->Fill( lInvMassXiMinus );
	if( lChargeXi > 0 && lIsBachelorMCPiPlus  )    fHistMassWithMcPIDXiPlus      ->Fill( lInvMassXiPlus );
	if( lChargeXi < 0 && lIsBachelorMCKMinus  )    fHistMassWithMcPIDOmegaMinus  ->Fill( lInvMassOmegaMinus );
	if( lChargeXi > 0 && lIsBachelorMCKPlus   )    fHistMassWithMcPIDOmegaPlus   ->Fill( lInvMassOmegaPlus );
	
        
	
	
	// - Step II.4 : MC association (care : lots of "continue;" below this line)
	//-------------	
	
	Bool_t lAssoXiMinus    = kFALSE;
	Bool_t lAssoXiPlus     = kFALSE;
	Bool_t lAssoOmegaMinus = kFALSE;
	Bool_t lAssoOmegaPlus  = kFALSE;
	
	
	if(fDebug > 5)
		cout 	<< "MC EventNumber : " << lMCevent->Header()->GetEvent() 
			<< " / MC event Number in Run : " << lMCevent->Header()->GetEventNrInRun() << endl;
	
	// - Step 4.1 : level of the V0 daughters
		
	Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );  
		// Abs value = needed ! question of quality track association ...
	Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );
		
	TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
	TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
	

	// - Step 4.2 : level of the Xi daughters
		
	Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ; 
	Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();
	
		if( lblMotherPosV0Dghter != lblMotherNegV0Dghter) continue; // same mother
		if( lblMotherPosV0Dghter < 0 ) continue; // mother != primary (!= -1)
		if( lblMotherNegV0Dghter < 0 ) continue;
					

		// mothers = Lambda candidate ... a priori
	
	TParticle* mcMotherPosV0Dghter = lMCstack->Particle( lblMotherPosV0Dghter );
	TParticle* mcMotherNegV0Dghter = lMCstack->Particle( lblMotherNegV0Dghter );

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
	
		
	// - Step 4.4 : Manage boolean for association
	
	if( mcMotherBach 		->GetPdgCode() ==   3312 &&
	    mcGdMotherPosV0Dghter	->GetPdgCode() ==   3312 &&
	    mcGdMotherNegV0Dghter	->GetPdgCode() ==   3312)	lAssoXiMinus = kTRUE;
	
	else if( mcMotherBach 		->GetPdgCode() ==  -3312 &&
	    mcGdMotherPosV0Dghter	->GetPdgCode() ==  -3312 &&
	    mcGdMotherNegV0Dghter	->GetPdgCode() ==  -3312)	lAssoXiPlus = kTRUE;
	
	else if( mcMotherBach 		->GetPdgCode() ==   3334 &&
	    mcGdMotherPosV0Dghter	->GetPdgCode() ==   3334 &&
	    mcGdMotherNegV0Dghter	->GetPdgCode() ==   3334)	lAssoOmegaMinus = kTRUE;
		
	else if( mcMotherBach 		->GetPdgCode() ==  -3334 &&
	    mcGdMotherPosV0Dghter	->GetPdgCode() ==  -3334 &&
	    mcGdMotherNegV0Dghter	->GetPdgCode() ==  -3334)	lAssoOmegaPlus = kTRUE;
	
	
	
	if(!lAssoXiMinus && !lAssoXiPlus && !lAssoOmegaMinus && !lAssoOmegaPlus) continue; // no association, skip the rest of the code	
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

	
	// - Step 5 : Plots around the cascade candidates associated with MC
	//-------------	
	
        Double_t lmcPt             = mcMotherBach->Pt();
        Double_t lmcRapCasc        = 0.5*TMath::Log( (mcMotherBach->Energy() + mcMotherBach->Pz()) / 
                                                     (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) );
        Double_t lmcPtBaryon       = (mcNegV0Dghter->GetPdgCode() == -2212 ) ? mcNegV0Dghter->Pt() : mcPosV0Dghter->Pt(); // spot the baryon daughter : pbar or p
        
        
        Double_t lmcEta            = mcMotherBach->Eta();
        Double_t lmcTransvRadius   = mcBach->R(); // to get the decay point of Xi, = the production vertex of Bachelor ...
        
        TVector3 lmcTVect3Mom( mcMotherBach->Px(), mcMotherBach->Py(), mcMotherBach->Pz() );

        Double_t lrecoPt           = xiESD->Pt();
        Double_t lrecoTransvRadius = TMath::Sqrt( xiESD->Xv() * xiESD->Xv() + xiESD->Yv() * xiESD->Yv() );
        
        TVector3 lrecoTVect3Mom( xiESD->Px(), xiESD->Py(), xiESD->Pz()  );
        Double_t lDeltaPhiMcReco   = lmcTVect3Mom.DeltaPhi( lrecoTVect3Mom ) * 180.0/TMath::Pi();

        
	// - Histos for the cascade candidates associated with MC
	
	if( lChargeXi < 0 && lAssoXiMinus){	
		fHistAsMCMassXiMinus	      ->Fill( lInvMassXiMinus  );
		if(lIsBachelorPion)	f2dHistAsMCandCombPIDGenPtVsGenYXiMinus->Fill( lmcPt, lmcRapCasc );
		f2dHistAsMCGenPtVsGenYXiMinus ->Fill( lmcPt, lmcRapCasc);
		fHistAsMCGenEtaXiMinus        ->Fill( lmcEta           );
		f2dHistAsMCResPtXiMinus       ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResRXiMinus        ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
                f2dHistAsMCResPhiXiMinus      ->Fill( lmcPt, lDeltaPhiMcReco );
                f2dHistAsMCPtProtonVsPtXiMinus->Fill( lmcPtBaryon, lmcPt );
	}
	
	else if( lChargeXi > 0 && lAssoXiPlus){	
                fHistAsMCMassXiPlus	                ->Fill( lInvMassXiPlus   );
                if(lIsBachelorPion)	f2dHistAsMCandCombPIDGenPtVsGenYXiPlus->Fill( lmcPt, lmcRapCasc );
                f2dHistAsMCGenPtVsGenYXiPlus            ->Fill( lmcPt, lmcRapCasc);
                fHistAsMCGenEtaXiPlus                   ->Fill( lmcEta           );
                f2dHistAsMCResPtXiPlus                  ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
                f2dHistAsMCResRXiPlus                   ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
                f2dHistAsMCResPhiXiPlus                 ->Fill( lmcPt, lDeltaPhiMcReco );
                f2dHistAsMCPtAntiProtonVsPtXiPlus       ->Fill( lmcPtBaryon, lmcPt );
    
	}
	
	else if( lChargeXi < 0 && lAssoOmegaMinus){	
		fHistAsMCMassOmegaMinus          ->Fill( lInvMassOmegaMinus );
		if(lIsBachelorKaon)	f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus->Fill( lmcPt, lmcRapCasc );
		f2dHistAsMCGenPtVsGenYOmegaMinus ->Fill( lmcPt, lmcRapCasc  );
		fHistAsMCGenEtaOmegaMinus        ->Fill( lmcEta             );
		f2dHistAsMCResPtOmegaMinus       ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResROmegaMinus        ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
                f2dHistAsMCResPhiOmegaMinus      ->Fill( lmcPt, lDeltaPhiMcReco );
                f2dHistAsMCPtProtonVsPtOmegaMinus->Fill( lmcPtBaryon, lmcPt );
	}
	
	else if( lChargeXi > 0 && lAssoOmegaPlus){	
                fHistAsMCMassOmegaPlus                  ->Fill( lInvMassOmegaPlus );
                if(lIsBachelorKaon)	f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus->Fill( lmcPt, lmcRapCasc );
                f2dHistAsMCGenPtVsGenYOmegaPlus         ->Fill( lmcPt, lmcRapCasc );
                fHistAsMCGenEtaOmegaPlus                ->Fill( lmcEta            );
                f2dHistAsMCResPtOmegaPlus               ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
                f2dHistAsMCResROmegaPlus                ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
                f2dHistAsMCResPhiOmegaPlus              ->Fill( lmcPt, lDeltaPhiMcReco );
                f2dHistAsMCPtAntiProtonVsPtOmegaPlus    ->Fill( lmcPtBaryon, lmcPt );
        }
        
        
        // - Step 6 : Containers = Cascade cuts + PID
	//-------------	

	// Double_t lChi2Xi         = -1. ;
        Double_t lDcaXiDaughters            = -1. ;
        Double_t lDcaBachToPrimVertexXi     = -1. ;
        Double_t lXiCosineOfPointingAngle   = -1. ;
        Double_t lPosXi[3]                  = { -1000.0, -1000.0, -1000.0 };
        Double_t lXiRadius                  = -1000. ;
        
        Double_t lInvMassLambdaAsCascDghter = 0.;
        Double_t lDcaV0DaughtersXi          = -1.;
        // Double_t lV0Chi2Xi               = -1. ;
        Double_t lV0CosineOfPointingAngleXi = -1.;
        Double_t lPosV0Xi[3]                = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
        Double_t lV0RadiusXi                = -1000.;
        Double_t lDcaV0ToPrimVertexXi       = -1.;
        Double_t lDcaPosToPrimVertexXi      = -1.;
        Double_t lDcaNegToPrimVertexXi      = -1.;

        //Int_t    nTrackWithTPCrefitMultiplicity  =  0;

        
        // 6.2 - Definition of the needed variables
        
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
        lXiRadius                       = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
        lInvMassLambdaAsCascDghter      = xiESD->GetEffMass(); 
                                        // This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
        lDcaV0DaughtersXi               = xiESD->GetDcaV0Daughters();
        // lV0Chi2Xi                    = xiESD->GetChi2V0();
        lV0CosineOfPointingAngleXi      = xiESD->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                             lBestPrimaryVtxPos[1],
                                                                             lBestPrimaryVtxPos[2] );
                xiESD->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] ); 
        lV0RadiusXi                     = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
        
        lDcaV0ToPrimVertexXi            = xiESD->GetD( lBestPrimaryVtxPos[0], 
                                                lBestPrimaryVtxPos[1], 
                                                lBestPrimaryVtxPos[2] );
        
        lDcaPosToPrimVertexXi           = TMath::Abs( pTrackXi   ->GetD( lBestPrimaryVtxPos[0],
                                                                         lBestPrimaryVtxPos[1],
                                                                         lMagneticField  )     );
        
        lDcaNegToPrimVertexXi           = TMath::Abs( nTrackXi   ->GetD( lBestPrimaryVtxPos[0], 
                                                                         lBestPrimaryVtxPos[1], 
                                                                         lMagneticField  )     );
        
        
        
        //nTrackWithTPCrefitMultiplicity  = DoESDTrackWithTPCrefitMultiplicity(lESDevent);  // FIXME : variable which is not used anymore at the moment ... 
                                                                                          //    ->  keep it while the task is still under development.
        
        
        
        // 6.3 - Filling the AliCFContainer (optimisation of topological selections + systematics)
        Double_t lContainerCutVars[20] = {0.0};

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

        lContainerCutVars[16] = lBestPrimaryVtxPos[2];
        lContainerCutVars[17] = nTrackPrimaryMultiplicity;       // FIXME : nTrackPrimaryMultiplicity not checked for AOD ...
        lContainerCutVars[18] = nITSandTPCtracksAndSPDtracklets; // FIXME : nITSandTPCtracksAndSPDtracklets is not available for AOD ... 
        lContainerCutVars[19] = lBachTPCClusters;                // FIXME : BachTPCClusters          is not available for AOD ... 

        // All cases should be covered below
        if( lChargeXi < 0 && lAssoXiMinus    ) {
                lContainerCutVars[11] = lInvMassXiMinus;
                lContainerCutVars[12] = 1.63;
                lContainerCutVars[14] = lmcRapCasc;
                lContainerCutVars[15] = -1.;
                        if( lIsBachelorPionForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )    
                                fCFContAsCascadeCuts->Fill(lContainerCutVars,0); // for Xi-
        }
        if( lChargeXi > 0 && lAssoXiPlus     ){
                lContainerCutVars[11] = lInvMassXiPlus;
                lContainerCutVars[12] = 1.26;
                lContainerCutVars[14] = lmcRapCasc;
                lContainerCutVars[15] = -1.; 
                        if( lIsBachelorPionForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )    
                                fCFContAsCascadeCuts->Fill(lContainerCutVars,1); // for Xi+
        }
        if( lChargeXi < 0 && lAssoOmegaMinus )  {
                lContainerCutVars[11] = 1.63;
                lContainerCutVars[12] = lInvMassOmegaMinus;
                lContainerCutVars[14] = -1.;
                lContainerCutVars[15] = lmcRapCasc;
                        if( lIsBachelorKaonForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC  && (TMath::Abs( lInvMassXiMinus-1.3217 ) > 0.008) )    
                                fCFContAsCascadeCuts->Fill(lContainerCutVars,2); // for Omega-
        }
	if( lChargeXi > 0 && lAssoOmegaPlus  ){
                lContainerCutVars[11] = 1.26;
                lContainerCutVars[12] = lInvMassOmegaPlus;
                lContainerCutVars[14] = -1.;
                lContainerCutVars[15] = lmcRapCasc;
                        if( lIsBachelorKaonForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC && (TMath::Abs( lInvMassXiPlus-1.3217 ) > 0.008) )    
                                fCFContAsCascadeCuts->Fill(lContainerCutVars,3); // for Omega+
        }
        
        
	// 6.4 - Filling the AliCFContainers related to PID

	Double_t lContainerPIDVars[4] = {0.0};

	
	// Xi Minus		
	if( lChargeXi < 0 && lAssoXiMinus ) {
		lContainerPIDVars[0] = lmcPt              ;
		lContainerPIDVars[1] = lInvMassXiMinus    ;
		lContainerPIDVars[2] = lmcRapCasc         ;
		lContainerPIDVars[3] = nTrackPrimaryMultiplicity ;  // FIXME : nTrackPrimaryMultiplicity is not checked for AOD ... 
			
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
		
		if(lIsBachelorPion     && 
		   lIsPosInXiProton && 
		   lIsNegInXiPion    )
		 	fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
	
	lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Xi Plus		
	if( lChargeXi > 0 && lAssoXiPlus ) {
		lContainerPIDVars[0] = lmcPt              ;
		lContainerPIDVars[1] = lInvMassXiPlus     ;
		lContainerPIDVars[2] = lmcRapCasc           ;
		lContainerPIDVars[3] = nTrackPrimaryMultiplicity ;    // FIXME : nTrackPrimaryMultiplicity is not checked for AOD ...  
			
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
		
		if(lIsBachelorPion     && 
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
		lContainerPIDVars[3] = nTrackPrimaryMultiplicity ;  // FIXME : nTrackPrimaryMultiplicity is not checked for AOD ... 
			
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
		lContainerPIDVars[3] = nTrackPrimaryMultiplicity ;  // FIXME : nTrackPrimaryMultiplicity is not checked for AOD ... 
			
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
		
		if(lIsBachelorKaon     && 
		   lIsNegInOmegaProton && 
		   lIsPosInOmegaPion    )
		 	fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
	}
        
	
	
}// End of loop over reconstructed cascades
 
 
 
 
  // Post output data.
 PostData(1, fListHistCascade);
}      



Int_t AliAnalysisTaskCheckPerformanceCascade::DoESDTrackWithTPCrefitMultiplicity(const AliESDEvent *lESDevent)
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
                    // Note : the event multiplicity = analysis on its own... See Jacek's or Jan Fiete's analysis on dN/d(pt) and dN/d(eta)

        }// end loop over all event tracks
        return  nTrackWithTPCrefitMultiplicity;
}







//________________________________________________________________________
void AliAnalysisTaskCheckPerformanceCascade::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
	
  TList *cRetrievedList = 0x0;
  cRetrievedList = (TList*)GetOutputData(1);
  if(!cRetrievedList){
	Printf("ERROR - AliAnalysisTaskCheckPerformanceCascade : ouput data container list not available\n");
	return;
  }	
	
  fHistMCTrackMultiplicity = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistMCTrackMultiplicity")  );
  if (!fHistMCTrackMultiplicity) {
    Printf("ERROR - AliAnalysisTaskCheckPerformanceCascade : fHistMCTrackMultiplicity not available");
    return;
  }
  
   
  TCanvas *canCheckPerformanceCascade = new TCanvas("AliAnalysisTaskCheckPerformanceCascade","Multiplicity",10,10,510,510);
  canCheckPerformanceCascade->cd(1)->SetLogy();

  fHistMCTrackMultiplicity->SetMarkerStyle(22);
  fHistMCTrackMultiplicity->DrawCopy("E");

}
