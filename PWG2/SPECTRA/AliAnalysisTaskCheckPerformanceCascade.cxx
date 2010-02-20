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
//            Modified : A.Maire Jan2010, antonin.maire@ires.in2p3.fr
//-----------------------------------------------------------------


#include <Riostream.h>

#include "TList.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TMath.h"

#include "AliLog.h"
#include "AliHeader.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliCFContainer.h"
#include "AliESDpid.h"
// #include "AliV0vertexer.h"
// #include "AliCascadeVertexer.h"
#include "AliESDEvent.h"
#include "AliESDcascade.h"

#include "AliAODEvent.h"

#include "AliAnalysisTaskCheckPerformanceCascade.h"

ClassImp(AliAnalysisTaskCheckPerformanceCascade)



     //_____Dummy constructor________________________________________________________________
AliAnalysisTaskCheckPerformanceCascade::AliAnalysisTaskCheckPerformanceCascade() 
: AliAnalysisTaskSE(), // <- take care to AliAnalysisTask( empty )
  fDebugCascade(0), fAnalysisType("ESD"), fCollidingSystems(0), fESDpid(0),
    
	// - Cascade part initialisation
fListHistCascade(0),
    fHistMCTrackMultiplicity(0), 

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

    fCFContCascadePIDAsXiMinus(0),
    fCFContCascadePIDAsXiPlus(0),
    fCFContCascadePIDAsOmegaMinus(0),
    fCFContCascadePIDAsOmegaPlus(0)

{
// Dummy constructor
}
     
       
     
     
//_____Non-default Constructor________________________________________________________________
AliAnalysisTaskCheckPerformanceCascade::AliAnalysisTaskCheckPerformanceCascade(const char *name) 
  : AliAnalysisTaskSE(name),
    fDebugCascade(0), fAnalysisType("ESD"), fCollidingSystems(0), fESDpid(0),
      
    	// - Cascade part initialisation
fListHistCascade(0),
    fHistMCTrackMultiplicity(0), 

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

    fCFContCascadePIDAsXiMinus(0),
    fCFContCascadePIDAsXiPlus(0),
    fCFContCascadePIDAsOmegaMinus(0),
    fCFContCascadePIDAsOmegaPlus(0)	

{
  // Constructor

  // Define input and output slots here
	// Input slot #0 works with a TChain
	//DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container (cascade)
  DefineOutput(1, TList::Class());
 
}





//________________________________________________________________________
void AliAnalysisTaskCheckPerformanceCascade::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

	
   // Option for AliLog
	AliLog::SetGlobalLogLevel(AliLog::kError); 
   	// to suppress the extensive info prompted by a run with MC			

if(! fESDpid){
		
  Double_t lAlephParameters[5] = {0.};
  	// Reasonable parameters extracted for p-p simulation (LHC09a4) - A.Kalweit
	// lAlephParameters[0] = 4.23232575531564326e+00/50;//50*0.76176e-1;
	// lAlephParameters[1] = 8.68482806165147636e+00;//10.632; 
	// lAlephParameters[2] = 1.34000000000000005e-05;//0.13279e-4;
	// lAlephParameters[3] = 2.30445734159456084e+00;//1.8631;
	// lAlephParameters[4] = 2.25624744086878559e+00;//1.9479;
        
        // Param for LHC09d10 prod - A.Kalweit
        lAlephParameters[0] = 2.15898e+00/50.;
        lAlephParameters[1] = 1.75295e+01;
        lAlephParameters[2] = 3.40030e-09;
        lAlephParameters[3] = 1.96178e+00;
        lAlephParameters[4] = 3.91720e+00; 
        
  fESDpid = new AliESDpid();
  fESDpid->GetTPCResponse().SetBetheBlochParameters(lAlephParameters[0],
					  lAlephParameters[1],
					  lAlephParameters[2],
					  lAlephParameters[3],
					  lAlephParameters[4]);
}
	
	
  // Definition of the datamembers	
  fListHistCascade = new TList();

  
  // - General
  
  if (!fHistMCTrackMultiplicity) {
     fHistMCTrackMultiplicity = new TH1F("fHistMCTrackMultiplicity", "MC Track Multiplicity;Number of MC tracks;Events", 100, 0, 500);
   //  fHistMCTrackMultiplicity = new TH1F("fHistMCTrackMultiplicity", "Multiplicity distribution;Number of tracks;Events", 200, 0, 40000); //HERE
  fListHistCascade->Add(fHistMCTrackMultiplicity);
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
     f2dHistGenPtVsGenYGenXiMinus = new TH2F("f2dHistGenPtVsGenYGenXiMinus", "MC P_{t} Vs MC Y of Gen #Xi^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
     fListHistCascade->Add(f2dHistGenPtVsGenYGenXiMinus);
  }
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascXiMinus) {
     fHistThetaGenCascXiMinus = new TH1F("fHistThetaGenCascXiMinus", "#theta of gen. #Xi^{-};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascXiMinus);
  }

  if (!f2dHistGenPtVsGenYFdblXiMinus) {
     f2dHistGenPtVsGenYFdblXiMinus = new TH2F("f2dHistGenPtVsGenYFdblXiMinus", "MC P_{t} Vs MC Y of findable Gen #Xi^{-}; Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
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
     fHistEtaGenCascXiPlus = new TH1F("fHistEtaGenCascXiPlus", "#eta of any gen. #Xi^{+};#eta;Number of Casc", 200, -10, 10);
     fListHistCascade->Add(fHistEtaGenCascXiPlus);
  }
  
  if (!f2dHistGenPtVsGenYGenXiPlus) {
     f2dHistGenPtVsGenYGenXiPlus = new TH2F("f2dHistGenPtVsGenYGenXiPlus", "MC P_{t} Vs MC Y of Gen #Xi^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
     fListHistCascade->Add(f2dHistGenPtVsGenYGenXiPlus);
  }
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascXiPlus) {
     fHistThetaGenCascXiPlus = new TH1F("fHistThetaGenCascXiPlus", "#theta of gen. #Xi^{+};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascXiPlus);
  }
 
  if (!f2dHistGenPtVsGenYFdblXiPlus) {
     f2dHistGenPtVsGenYFdblXiPlus = new TH2F("f2dHistGenPtVsGenYFdblXiPlus", "MC P_{t} Vs MC Y of findable Gen #Xi^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
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
     f2dHistGenPtVsGenYGenOmegaMinus = new TH2F("f2dHistGenPtVsGenYGenOmegaMinus", "MC P_{t} Vs MC Y of Gen #Omega^{-} ;Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
     fListHistCascade->Add(f2dHistGenPtVsGenYGenOmegaMinus);
  }
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascOmegaMinus) {
     fHistThetaGenCascOmegaMinus = new TH1F("fHistThetaGenCascOmegaMinus", "#theta of gen. #Omega^{-};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascOmegaMinus);
  }
 
  if (!f2dHistGenPtVsGenYFdblOmegaMinus) {
     f2dHistGenPtVsGenYFdblOmegaMinus = new TH2F("f2dHistGenPtVsGenYFdblOmegaMinus", "MC P_{t} Vs MC Y of findable Gen #Omega^{-}; Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
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
  
  if (!f2dHistGenPtVsGenYGenOmegaPlus) {
     f2dHistGenPtVsGenYGenOmegaPlus = new TH2F("f2dHistGenPtVsGenYGenOmegaPlus", "MC P_{t} Vs MC Y of Gen #Omega^{+} ;Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
     fListHistCascade->Add(f2dHistGenPtVsGenYGenOmegaPlus);
  }
  
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascOmegaPlus) {
     fHistThetaGenCascOmegaPlus = new TH1F("fHistThetaGenCascOmegaPlus", "#theta of gen. #Omega^{+};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascOmegaPlus);
  }
 
  if (!f2dHistGenPtVsGenYFdblOmegaPlus) {
     f2dHistGenPtVsGenYFdblOmegaPlus = new TH2F("f2dHistGenPtVsGenYFdblOmegaPlus", "MC P_{t} Vs MC Y of findable Gen #Omega^{+}; Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
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
	  fHistMassXiMinus = new TH1F("fHistMassXiMinus","#Xi^{-} candidates;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 200,1.2,2.0);
	  fListHistCascade->Add(fHistMassXiMinus);
  }
  
  if (! fHistMassXiPlus) {
	  fHistMassXiPlus = new TH1F("fHistMassXiPlus","#Xi^{+} candidates;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",200,1.2,2.0);
	  fListHistCascade->Add(fHistMassXiPlus);
  }

  if (! fHistMassOmegaMinus) {
	  fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus","#Omega^{-} candidates;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 250,1.5,2.5);
	  fListHistCascade->Add(fHistMassOmegaMinus);
  }
 
  if (! fHistMassOmegaPlus) {
	  fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus","#Omega^{+} candidates;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 250,1.5,2.5);
	  fListHistCascade->Add(fHistMassOmegaPlus);
  }
  
  
  
  		// - Effective mass histos with combined PID
  
  if (! fHistMassWithCombPIDXiMinus) {
    fHistMassWithCombPIDXiMinus = new TH1F("fHistMassWithCombPIDXiMinus","#Xi^{-} candidates, with Bach. comb. PID;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 200,1.2,2.0);
    fListHistCascade->Add(fHistMassWithCombPIDXiMinus);
  }
  
  if (! fHistMassWithCombPIDXiPlus) {
    fHistMassWithCombPIDXiPlus = new TH1F("fHistMassWithCombPIDXiPlus","#Xi^{+} candidates, with Bach. comb. PID;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",200,1.2,2.0);
    fListHistCascade->Add(fHistMassWithCombPIDXiPlus);
  }

  if (! fHistMassWithCombPIDOmegaMinus) {
	fHistMassWithCombPIDOmegaMinus = new TH1F("fHistMassWithCombPIDOmegaMinus","#Omega^{-} candidates, with Bach. comb. PID;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 250,1.5,2.5);
    fListHistCascade->Add(fHistMassWithCombPIDOmegaMinus);
  }
 
  if (! fHistMassWithCombPIDOmegaPlus) {
	fHistMassWithCombPIDOmegaPlus = new TH1F("fHistMassWithCombPIDOmegaPlus","#Omega^{+} candidates, with Bach. comb. PID;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 250,1.5,2.5);
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
    fHistMassWithMcPIDXiMinus = new TH1F("fHistMassWithMcPIDXiMinus","#Xi^{-} candidates, with Bach. MC PID;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 200,1.2,2.0);
    fListHistCascade->Add(fHistMassWithMcPIDXiMinus);
  }
  
  if (! fHistMassWithMcPIDXiPlus) {
    fHistMassWithMcPIDXiPlus = new TH1F("fHistMassWithMcPIDXiPlus","#Xi^{+} candidates, with Bach. MC PID;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",200,1.2,2.0);
    fListHistCascade->Add(fHistMassWithMcPIDXiPlus);
  }

  if (! fHistMassWithMcPIDOmegaMinus) {
	fHistMassWithMcPIDOmegaMinus = new TH1F("fHistMassWithMcPIDOmegaMinus","#Omega^{-} candidates, with Bach. MC PID;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 250,1.5,2.5);
    fListHistCascade->Add(fHistMassWithMcPIDOmegaMinus);
  }
 
  if (! fHistMassWithMcPIDOmegaPlus) {
	fHistMassWithMcPIDOmegaPlus = new TH1F("fHistMassWithMcPIDOmegaPlus","#Omega^{+} candidates, with Bach. MC PID;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 250,1.5,2.5);
    fListHistCascade->Add(fHistMassWithMcPIDOmegaPlus);
  }
  
  
		// - Effective mass histos for cascades candidates ASSOCIATED with MC.
  
  if (! fHistAsMCMassXiMinus) {
	  fHistAsMCMassXiMinus = new TH1F("fHistAsMCMassXiMinus","#Xi^{-} candidates associated to MC;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 200,1.2,2.0);
	  fListHistCascade->Add(fHistAsMCMassXiMinus);
  }
  
  if (! fHistAsMCMassXiPlus) {
	  fHistAsMCMassXiPlus = new TH1F("fHistAsMCMassXiPlus","#Xi^{+} candidates associated to MC;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",200,1.2,2.0);
	  fListHistCascade->Add(fHistAsMCMassXiPlus);
  }

  if (! fHistAsMCMassOmegaMinus) {
	  fHistAsMCMassOmegaMinus = new TH1F("fHistAsMCMassOmegaMinus","#Omega^{-} candidates associated to MC;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 250,1.5,2.5);
	  fListHistCascade->Add(fHistAsMCMassOmegaMinus);
  }
 
  if (! fHistAsMCMassOmegaPlus) {
	  fHistAsMCMassOmegaPlus = new TH1F("fHistAsMCMassOmegaPlus","#Omega^{+} candidates associated to MC;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 250,1.5,2.5);
	  fListHistCascade->Add(fHistAsMCMassOmegaPlus);
  }
  
		
		// -  Generated Pt Vs generated Y of the cascade candidates associated with MC 
		//     + having the proper maximum proba of combined PID for the bachelor
  
  if (!f2dHistAsMCandCombPIDGenPtVsGenYXiMinus) {
     f2dHistAsMCandCombPIDGenPtVsGenYXiMinus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYXiMinus", "MC P_{t} Vs MC Y of #Xi^{-} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
     fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYXiMinus);
  }
  
  if (!f2dHistAsMCandCombPIDGenPtVsGenYXiPlus) {
     f2dHistAsMCandCombPIDGenPtVsGenYXiPlus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYXiPlus", "MC P_{t} Vs MC Y of #Xi^{+} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
     fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYXiPlus);
  }
  
  if (!f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus) {
     f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus", "MC P_{t} Vs MC Y of #Omega^{-} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
     fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus);
  }
  
  if (!f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus) {
     f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus = new TH2F("f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus", "MC P_{t} Vs MC Y of #Omega^{+} (associated+Bach.PID); Pt_{MC} (GeV/c); Y_{MC}", 200, 0., 10., 240, -1.2, 1.2);
     fListHistCascade->Add(f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus);
  }
  

  		// - Generated Pt Vs Generated Y, for the cascade candidates associated with MC
  
  if (!f2dHistAsMCGenPtVsGenYXiMinus) {
  	f2dHistAsMCGenPtVsGenYXiMinus = new TH2F("f2dHistAsMCGenPtVsGenYXiMinus", "MC P_{t} Vs MC Y of gen. #Xi^{-} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 240, -1.2, 1.2);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYXiMinus );
  }
  
  if (!f2dHistAsMCGenPtVsGenYXiPlus) {
	  f2dHistAsMCGenPtVsGenYXiPlus = new TH2F("f2dHistAsMCGenPtVsGenYXiPlus", "MC P_{t} Vs MC Y of gen. #Xi^{+} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 240, -1.2, 1.2);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYXiPlus );
  }
  
  if (!f2dHistAsMCGenPtVsGenYOmegaMinus) {
	  f2dHistAsMCGenPtVsGenYOmegaMinus = new TH2F("f2dHistAsMCGenPtVsGenYOmegaMinus", "MC P_{t} Vs MC Y of gen. #Omega^{-} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 240, -1.2, 1.2);
	  fListHistCascade->Add(f2dHistAsMCGenPtVsGenYOmegaMinus );
  }
  
  if (!f2dHistAsMCGenPtVsGenYOmegaPlus) {
	  f2dHistAsMCGenPtVsGenYOmegaPlus = new TH2F("f2dHistAsMCGenPtVsGenYOmegaPlus", "MC P_{t} Vs MC Y of gen. #Omega^{+} (associated);Pt_{MC} (GeV/c); Rapidity, Y_{MC}",200, 0., 10., 240, -1.2, 1.2);
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
  
  
  
  
                // - PID container
if(! fCFContCascadePIDAsXiMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4];
  lNbBinsPerVar[0] = 200;
  lNbBinsPerVar[1] = 400;
  lNbBinsPerVar[2] = 48;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDAsXiMinus = new AliCFContainer("fCFContCascadePIDAsXiMinus","Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDAsXiMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDAsXiMinus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
  fCFContCascadePIDAsXiMinus->SetBinLimits(2,  -1.2  ,   1.2 );	// Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDAsXiMinus->SetBinLimits(3, 0.0, 20000.0  );    // TPCrefitTrackMultiplicity
  else
	fCFContCascadePIDAsXiMinus->SetBinLimits(3, 0.0, 250.0  );     // TPCrefitTrackMultiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDAsXiMinus->SetStepTitle(0, "No PID");
  fCFContCascadePIDAsXiMinus->SetStepTitle(1, "TPC PID / 3-#sigma cut on Bachelor track");
  fCFContCascadePIDAsXiMinus->SetStepTitle(2, "TPC PID / 3-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDAsXiMinus->SetStepTitle(3, "TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDAsXiMinus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDAsXiMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDAsXiMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDAsXiMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDAsXiMinus->SetVarTitle(1, "M( #Lambda , #pi^{-} ) (GeV/c^{2})");
  fCFContCascadePIDAsXiMinus->SetVarTitle(2, "Y_{#Xi}");
  fCFContCascadePIDAsXiMinus->SetVarTitle(3, "TPCrefit track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDAsXiMinus);
  
}

if(! fCFContCascadePIDAsXiPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4];
  lNbBinsPerVar[0] = 200;
  lNbBinsPerVar[1] = 400;
  lNbBinsPerVar[2] = 48;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDAsXiPlus = new AliCFContainer("fCFContCascadePIDAsXiPlus","Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDAsXiPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDAsXiPlus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
  fCFContCascadePIDAsXiPlus->SetBinLimits(2,  -1.2  ,   1.2 );	// Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDAsXiPlus->SetBinLimits(3, 0.0, 20000.0  );    // TPCrefitTrackMultiplicity
  else
	fCFContCascadePIDAsXiPlus->SetBinLimits(3, 0.0, 250.0  );     // TPCrefitTrackMultiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDAsXiPlus->SetStepTitle(0, "No PID");
  fCFContCascadePIDAsXiPlus->SetStepTitle(1, "TPC PID / 3-#sigma cut on Bachelor track");
  fCFContCascadePIDAsXiPlus->SetStepTitle(2, "TPC PID / 3-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDAsXiPlus->SetStepTitle(3, "TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDAsXiPlus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDAsXiPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDAsXiPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDAsXiPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDAsXiPlus->SetVarTitle(1, "M( #Lambda , #pi^{+} ) (GeV/c^{2})");
  fCFContCascadePIDAsXiPlus->SetVarTitle(2, "Y_{#Xi}");
  fCFContCascadePIDAsXiPlus->SetVarTitle(3, "TPCrefit track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDAsXiPlus);
  
}


if(! fCFContCascadePIDAsOmegaMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4];
  lNbBinsPerVar[0] = 200;
  lNbBinsPerVar[1] = 500;
  lNbBinsPerVar[2] = 48;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDAsOmegaMinus = new AliCFContainer("fCFContCascadePIDAsOmegaMinus","Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDAsOmegaMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDAsOmegaMinus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
  fCFContCascadePIDAsOmegaMinus->SetBinLimits(2,  -1.2  ,   1.2 );	// Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDAsOmegaMinus->SetBinLimits(3, 0.0, 20000.0  );    // TPCrefitTrackMultiplicity
  else
	fCFContCascadePIDAsOmegaMinus->SetBinLimits(3, 0.0, 250.0  );     // TPCrefitTrackMultiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(0, "No PID");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(1, "TPC PID / 3-#sigma cut on Bachelor track");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(2, "TPC PID / 3-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(3, "TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDAsOmegaMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDAsOmegaMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDAsOmegaMinus->SetVarTitle(1, "M( #Lambda , K^{-} ) (GeV/c^{2})");
  fCFContCascadePIDAsOmegaMinus->SetVarTitle(2, "Y_{#Omega}");
  fCFContCascadePIDAsOmegaMinus->SetVarTitle(3, "TPCrefit track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDAsOmegaMinus);
  
}

if(! fCFContCascadePIDAsOmegaPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4];
  lNbBinsPerVar[0] = 200;
  lNbBinsPerVar[1] = 500;
  lNbBinsPerVar[2] = 48;
  lNbBinsPerVar[3] = 250;
   
  
  fCFContCascadePIDAsOmegaPlus = new AliCFContainer("fCFContCascadePIDAsOmegaPlus","Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDAsOmegaPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDAsOmegaPlus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
  fCFContCascadePIDAsOmegaPlus->SetBinLimits(2,  -1.2  ,   1.2 );	// Rapidity
  if(fCollidingSystems) 
	fCFContCascadePIDAsOmegaPlus->SetBinLimits(3, 0.0, 20000.0  );    // TPCrefitTrackMultiplicity
  else
	fCFContCascadePIDAsOmegaPlus->SetBinLimits(3, 0.0, 250.0  );     // TPCrefitTrackMultiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(0, "No PID");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(1, "TPC PID / 3-#sigma cut on Bachelor track");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(2, "TPC PID / 3-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(3, "TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDAsOmegaPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDAsOmegaPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDAsOmegaPlus->SetVarTitle(1, "M( #Lambda , K^{+} ) (GeV/c^{2})");
  fCFContCascadePIDAsOmegaPlus->SetVarTitle(2, "Y_{#Omega}");
  fCFContCascadePIDAsOmegaPlus->SetVarTitle(3, "TPCrefit track Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDAsOmegaPlus);
  
}

  
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
			cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
			return;
		}
		ncascades = lESDevent->GetNumberOfCascades();
	}
  
	else if(fAnalysisType == "AOD"){  
		lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() ); 
		if (!lAODevent) {
			Printf("ERROR: lAODevent not available \n");
			cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
			return;
		}
		ncascades = lAODevent->GetNumberOfCascades();
	}
	

	lMCevent = MCEvent();
	if (!lMCevent) {
		Printf("ERROR: Could not retrieve MC event \n");
		cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;	
		return;
	}

	lMCstack = lMCevent->Stack();
	if (!lMCstack) {
		Printf("ERROR: Could not retrieve MC stack \n");
		cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
		return;
		
	}
	
  //-------------------------------------------------
  // 0 - Trigger managment
  // FIXME : Check the availability of the proper trigger 

   // Note : Presuppose the presence of AliPhysicsSelectionTask
        Bool_t isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
        if ( ! isSelected ) return;
        //else Printf("Event selected ... \n");
        
        
  //	cout << "Name of the accessed file :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;

  //	cout << "Tree characteristics ..." << endl;
  //	fInputHandler->GetTree()->Print("toponly");
  //	fInputHandler->GetTree()->GetBranch("PrimaryVertex")->Print();
  //	fInputHandler->GetTree()->GetBranch("SPDVertex")->Print();
        
        
   //-------------------------------------------------
   // 1 - Cascade vertexer (ESD)

        // if(fAnalysisType == "ESD" ){
        //         lESDevent->ResetCascades();
        //         lESDevent->ResetV0s();
        // 
        //         AliV0vertexer V0vtxer;
        //         AliCascadeVertexer CascVtxer;
        // 
        //         Double_t v0sels[]={33,    // max allowed chi2
        //                         0.02,  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
        //                         0.02,  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
        //                         1.0,   // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
        //                         0.98,  // max allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
        //                         0.2,   // min radius of the fiducial volume                 (LHC09a4 : 0.2)
        //                         100    // max radius of the fiducial volume                 (LHC09a4 : 100.0)
        //         };
        //         V0vtxer.SetDefaultCuts(v0sels);
        // 
        //         Double_t xisels[]={33.,   // max allowed chi2 (same as PDC07)
        //                         0.02,   // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
        //                         0.020,  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
        //                         0.01 ,  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
        //                         0.5,    // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
        //                         0.99,   // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
        //                         0.2,    // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
        //                         100     // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
        //         };
        //         CascVtxer.SetDefaultCuts(xisels);
        // 
        //         V0vtxer.Tracks2V0vertices(lESDevent);
        //         CascVtxer.V0sTracks2CascadeVertices(lESDevent);
        // }


 

 
  // ---------------------------------------------------------------
  // - Initialisation of the part dedicated to cascade vertices
  
  Int_t iNumberOfPrimaries = -1;
  iNumberOfPrimaries = lMCstack->GetNprimary();
  
	if(iNumberOfPrimaries < 1) return;
    
       fHistMCTrackMultiplicity->Fill( lMCstack->GetNtrack() );
   
   // For proton
   /*
   for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < iNumberOfPrimaries; iCurrentLabelStack++) 
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

      
       
//__________________________________________________________________________	
// Part 1 - Loop over the different types of generated cascades (Xi-+, Omega-+)	

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


   for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < iNumberOfPrimaries; iCurrentLabelStack++) 
    {// This is the begining of the loop on primaries
      
        TParticle* lCurrentParticle = 0x0; 
    	           lCurrentParticle = lMCstack->Particle( iCurrentLabelStack );
	if(!lCurrentParticle){
		Printf("Cascade loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
		continue;
		
	}
	 
    	if( lCurrentParticle->GetPdgCode() == lPdgCodeCasc ){  // Here !
   		//cout << "Xi- within loop " << iCurrentLabelStack << "/ " << iNumberOfPrimaries << endl;
		
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
		
			if( lBach->Pt() < 0.2 ) continue; //FIXME : maybe tuned for Xi but not for K- from Omega ...
			
		
		
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
			
			if( lDghtBarV0->Pt() < 0.3 ) continue;
			if( lDghtMesV0->Pt() < 0.2 ) continue;
			
			
			
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
// Part 2 - Loop over the reconstructed candidates
  
  
// Temporary way : AOD awareness of the code to be developed  
if(fAnalysisType == "AOD") return;


for (Int_t iXi = 0; iXi < ncascades; iXi++) 
{// This is the begining of the Cascade loop
		
	AliESDcascade *xiESD = lESDevent->GetCascade(iXi);
	if (!xiESD) continue;
	
	// -  Step 1 : Preparing the general info about of the event = prim. Vtx + magnetic field (ESD)
	//-------------
		
	// Double_t lMagneticField = lESDevent->GetMagneticField( );

        const AliESDVertex *lPrimaryTrackingVtx = lESDevent->GetPrimaryVertexTracks();  // get the vtx stored in ESD found with tracks
        const AliESDVertex *lPrimarySPDVtx      = lESDevent->GetPrimaryVertexSPD();     // get the vtx stored in ESD found with SPD tracklets

        const AliESDVertex *lPrimaryBestVtx     = lESDevent->GetPrimaryVertex();
                // get the best primary vertex available for the event
                // As done in AliCascadeVertexer, we keep the one which is the best one available.
                // between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
                Double_t lBestPrimaryVtxPos[3]   = {-100.0, -100.0, -100.0};
        lPrimaryBestVtx->GetXYZ( lBestPrimaryVtxPos );

        // FIXME : quality cut on the z-position of the prim vertex.
        Bool_t kQualityCutZprimVtxPos = kTRUE;
        if(kQualityCutZprimVtxPos) {
                if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0 ) { AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !"); return; }
        }
        // FIXME : remove TPC-only primary vertex : retain only events with tracking + SPD vertex
        Bool_t kQualityCutNoTPConlyPrimVtx = kTRUE;
        if(kQualityCutNoTPConlyPrimVtx) {
                if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingVtx->GetStatus() ){
                        AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
                        return;
                }
        }
	
	
	// - Step 2 : Connection to daughter tracks of the current cascade
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
        Bool_t kQualityCutTPCrefit = kFALSE;
        Bool_t kQualityCut80TPCcls = kFALSE;
        
        if(kQualityCutTPCrefit){
                // 1 - Poor quality related to TPCrefit
                ULong_t pStatus    = pTrackXi->GetStatus();
                ULong_t nStatus    = nTrackXi->GetStatus();
                ULong_t bachStatus = bachTrackXi->GetStatus();
                if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
                if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
                if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
        }
        if(kQualityCut80TPCcls){
                // 2 - Poor quality related to TPC clusters
                if(lPosTPCClusters  < 80) { AliWarning("Pb / V0 Pos. track has less than 80 TPC clusters ... continue!"); continue; }
                if(lNegTPCClusters  < 80) { AliWarning("Pb / V0 Neg. track has less than 80 TPC clusters ... continue!"); continue; }
                if(lBachTPCClusters < 80) { AliWarning("Pb / Bach.   track has less than 80 TPC clusters ... continue!"); continue; }
        }
	
	// - Step 3 : Info over reconstructed cascades
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
	
	
	// - Step 4 : PID info
	//-------------
	
	
	// 4.1 - PID Information

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

        // 4.1.A - Combined PID
	// Reasonable guess for the priors for the cascade track sample (e-, mu, pi, K, p)
	Double_t lPriorsGuessXi[5]    = {0, 0, 2, 0, 1};
	Double_t lPriorsGuessOmega[5] = {0, 0, 1, 1, 1};
	
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
	AliPID bachPidXi;	bachPidXi.SetPriors(    lPriorsGuessXi    );
	AliPID bachPidOmega;	bachPidOmega.SetPriors( lPriorsGuessOmega );
	
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
	
	
	// 4.1.B - TPC PID : 3-sigma bands on Bethe-Bloch curve
	// Bachelor
	if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 3) lIsBachelorKaonForTPC = kTRUE;
	if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 3) lIsBachelorPionForTPC = kTRUE;
	
	// Negative V0 daughter
	if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 3) lIsNegPionForTPC   = kTRUE;
	if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 3) lIsNegProtonForTPC = kTRUE;
	
	// Positive V0 daughter
	if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 3) lIsPosPionForTPC   = kTRUE;
	if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 3) lIsPosProtonForTPC = kTRUE;

		
        // Combined PID TH1s
	if( lChargeXi < 0 && lIsBachelorPion )    fHistMassWithCombPIDXiMinus     ->Fill( lInvMassXiMinus    );
	if( lChargeXi > 0 && lIsBachelorPion )    fHistMassWithCombPIDXiPlus      ->Fill( lInvMassXiPlus     );
	if( lChargeXi < 0 && lIsBachelorKaon )    fHistMassWithCombPIDOmegaMinus  ->Fill( lInvMassOmegaMinus );
	if( lChargeXi > 0 && lIsBachelorKaon )    fHistMassWithCombPIDOmegaPlus   ->Fill( lInvMassOmegaPlus  );
         
	
	// 4.2 - PID proba Vs Pt(Bach)
	Int_t      lblBachForPID  = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
	TParticle* mcBachForPID   = lMCstack->Particle( lblBachForPID );
	Double_t   lmcPtBach      = mcBachForPID->Pt();
	
	if(lIsBachelorPion)   f2dHistPIDprobaPionVsMCPtBach->Fill( lmcPtBach, ppionBach );
	if(lIsBachelorKaon)   f2dHistPIDprobaKaonVsMCPtBach->Fill( lmcPtBach, pkaonBach );
	
			
	// 4.3 - MC perfect PID
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
	
        
	
	
	// - Step 5 : MC association (care : lots of "continue;" below this line)
	//-------------	
	
	Bool_t lAssoXiMinus    = kFALSE;
	Bool_t lAssoXiPlus     = kFALSE;
	Bool_t lAssoOmegaMinus = kFALSE;
	Bool_t lAssoOmegaPlus  = kFALSE;
	
	
	if(fDebug > 5)
		cout 	<< "MC EventNumber : " << lMCevent->Header()->GetEvent() 
			<< " / MC event Number in Run : " << lMCevent->Header()->GetEventNrInRun() << endl;
	
	// - Step 5.1 : level of the V0 daughters
		
	Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );  
		// Abs value = needed ! question of quality track association ...
	Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );
		
	TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
	TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
	

	// - Step 5.2 : level of the Xi daughters
		
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
				

	// - Step 5.3 : level of Xi candidate
	
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
	
		
	// - Step 5.4 : Manage boolean for association
	
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
	
	
	
	if(!lAssoXiMinus && !lAssoXiPlus && !lAssoOmegaMinus && !lAssoOmegaPlus) continue; // no association
	
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

	
	// - Step 6 : Plots around the cascade candidates associated with MC
	//-------------	
	
        Double_t lmcPt             = mcMotherBach->Pt();
        Double_t lmcRapXi          = 0.5*TMath::Log( (mcMotherBach->Energy() + mcMotherBach->Pz()) / 
                                                     (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) );
        Double_t lmcEta            = mcMotherBach->Eta();
        Double_t lmcTransvRadius   = mcBach->R(); // to get the decay point of Xi, = the production vertex of Bachelor ...

        Double_t lrecoPt           = xiESD->Pt();
        Double_t lrecoTransvRadius = TMath::Sqrt( xiESD->Xv() * xiESD->Xv() + xiESD->Yv() * xiESD->Yv() );
	
	// - Histos for the cascade candidates associated with MC
	
	if( lChargeXi < 0 && lAssoXiMinus){	
		fHistAsMCMassXiMinus	      ->Fill( lInvMassXiMinus  );
		if(lIsBachelorPion)	f2dHistAsMCandCombPIDGenPtVsGenYXiMinus->Fill( lmcPt, lmcRapXi );
		f2dHistAsMCGenPtVsGenYXiMinus ->Fill( lmcPt, lmcRapXi  );
		fHistAsMCGenEtaXiMinus        ->Fill( lmcEta           );
		f2dHistAsMCResPtXiMinus       ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResRXiMinus        ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
	}
	
	else if( lChargeXi > 0 && lAssoXiPlus){	
		fHistAsMCMassXiPlus	      ->Fill( lInvMassXiPlus   );
		if(lIsBachelorPion)	f2dHistAsMCandCombPIDGenPtVsGenYXiPlus->Fill( lmcPt, lmcRapXi );
		f2dHistAsMCGenPtVsGenYXiPlus  ->Fill( lmcPt, lmcRapXi  );
		fHistAsMCGenEtaXiPlus         ->Fill( lmcEta           );
		f2dHistAsMCResPtXiPlus        ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResRXiPlus         ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
	}
	
	else if( lChargeXi < 0 && lAssoOmegaMinus){	
		fHistAsMCMassOmegaMinus          ->Fill( lInvMassOmegaMinus );
		if(lIsBachelorKaon)	f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus->Fill( lmcPt, lmcRapXi );
		f2dHistAsMCGenPtVsGenYOmegaMinus ->Fill( lmcPt, lmcRapXi    );
		fHistAsMCGenEtaOmegaMinus        ->Fill( lmcEta             );
		f2dHistAsMCResPtOmegaMinus       ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResROmegaMinus        ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
	}
	
	else if( lChargeXi > 0 && lAssoOmegaPlus){	
		fHistAsMCMassOmegaPlus           ->Fill( lInvMassOmegaPlus );
		if(lIsBachelorKaon)	f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus->Fill( lmcPt, lmcRapXi );
		f2dHistAsMCGenPtVsGenYOmegaPlus  ->Fill( lmcPt, lmcRapXi   );
		fHistAsMCGenEtaOmegaPlus         ->Fill( lmcEta            );
		f2dHistAsMCResPtOmegaPlus        ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResROmegaPlus         ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
	}
        
        
        
        
        
        //FIXME
        
	// - Filling the AliCFContainers related to PID
	
	Double_t lContainerPIDVars[4] = {0.0};
        Int_t nTrackWithTPCrefitMultiplicity =  0;
        
        nTrackWithTPCrefitMultiplicity = DoESDTrackWithTPCrefitMultiplicity(lESDevent);
        
        
        
        
	
	// Xi Minus		
	if( lChargeXi < 0 && lAssoXiMinus ) {
		lContainerPIDVars[0] = lmcPt              ;
		lContainerPIDVars[1] = lInvMassXiMinus    ;
		lContainerPIDVars[2] = lmcRapXi           ;
		lContainerPIDVars[3] = nTrackWithTPCrefitMultiplicity ;
			
		// No PID
			fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorPionForTPC  )
			fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 1); // TPC PID / 3-#sigma cut on Bachelor track
		
		if( lIsBachelorPionForTPC && 
		    lIsPosProtonForTPC     )
			fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 2); // TPC PID / 3-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorPionForTPC && 
		    lIsPosProtonForTPC    && 
		    lIsNegPionForTPC       )
			fCFContCascadePIDAsXiMinus->Fill(lContainerPIDVars, 3); // TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks
		
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
		lContainerPIDVars[2] = lmcRapXi           ;
		lContainerPIDVars[3] = nTrackWithTPCrefitMultiplicity ;
			
		// No PID
			fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorPionForTPC  )
			fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 1); // TPC PID / 3-#sigma cut on Bachelor track
		
		if( lIsBachelorPionForTPC && 
		    lIsNegProtonForTPC     )
			fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 2); // TPC PID / 3-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorPionForTPC && 
		    lIsNegProtonForTPC    && 
		    lIsPosPionForTPC       )
			fCFContCascadePIDAsXiPlus->Fill(lContainerPIDVars, 3); // TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks
		
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
		lContainerPIDVars[2] = lmcRapXi           ;
		lContainerPIDVars[3] = nTrackWithTPCrefitMultiplicity ;
			
		// No PID
			fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorKaonForTPC  )
			fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 1); // TPC PID / 3-#sigma cut on Bachelor track
		
		if( lIsBachelorKaonForTPC && 
		    lIsPosProtonForTPC     )
			fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 2); // TPC PID / 3-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorKaonForTPC && 
		    lIsPosProtonForTPC    && 
		    lIsNegPionForTPC       )
			fCFContCascadePIDAsOmegaMinus->Fill(lContainerPIDVars, 3); // TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks
		
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
		lContainerPIDVars[2] = lmcRapXi           ;
		lContainerPIDVars[3] = nTrackWithTPCrefitMultiplicity ;
			
		// No PID
			fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
		if( lIsBachelorKaonForTPC  )
			fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 1); // TPC PID / 3-#sigma cut on Bachelor track
		
		if( lIsBachelorKaonForTPC && 
		    lIsNegProtonForTPC     )
			fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 2); // TPC PID / 3-#sigma cut on Bachelor+Baryon tracks
		
		if( lIsBachelorKaonForTPC && 
		    lIsNegProtonForTPC    && 
		    lIsPosPionForTPC       )
			fCFContCascadePIDAsOmegaPlus->Fill(lContainerPIDVars, 3); // TPC PID / 3-#sigma cut on Bachelor+Baryon+Meson tracks
		
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
                    // Note : the event multiplicity = analysis on its own... See Jacek's or Jan Fiete's analysis on dN/d(eta)

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
