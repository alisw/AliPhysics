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
//                 AliAnalysisTaskCheckPerformanceCascade class
//            This task is for a performance study of cascade identification.
//            It works with MC info and ESD/AOD tree.
//            Origin : A.Maire Mar2009, antonin.maire@ires.in2p3.fr
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

#include "AliESDEvent.h"
#include "AliESDcascade.h"

#include "AliAODEvent.h"

#include "AliAnalysisTaskCheckPerformanceCascade.h"

ClassImp(AliAnalysisTaskCheckPerformanceCascade)



     //_____Dummy constructor________________________________________________________________
AliAnalysisTaskCheckPerformanceCascade::AliAnalysisTaskCheckPerformanceCascade() 
: AliAnalysisTaskSE(), // <- take care to AliAnalysisTask( empty )
  fDebugCascade(0), fAnalysisType("ESD"), fCollidingSystems(0),
    
	// - Cascade part initialisation
fListHistCascade(0),
    fHistMCTrackMultiplicity(0), 

    fHistEtaGenProton(0),
    fHistEtaGenAntiProton(0),
    
   // Xi-
   fHistEtaGenCascXiMinus(0),
      
   fHistYGenCascMidRapXiMinus(0),		
   fHistEtaGenCascMidRapXiMinus(0),
   fHistThetaGenCascMidRapXiMinus(0),
   fHistPtGenCascMidRapXiMinus(0),
   
    fHistThetaGenCascXiMinus(0),
    fHistPtFdblGenCascXiMinus(0),
    
    fHistThetaLambdaXiMinus(0), 
    fHistThetaBachXiMinus(0),
    
    fHistThetaMesDghterXiMinus(0), 
    fHistThetaBarDghterXiMinus(0),
    
    fHistPtBachXiMinus(0),
    fHistPtMesDghterXiMinus(0),
    fHistPtBarDghterXiMinus(0),
   
   
   // Xi+
   fHistEtaGenCascXiPlus(0),
      
    fHistYGenCascMidRapXiPlus(0),
    fHistEtaGenCascMidRapXiPlus(0),
    fHistThetaGenCascMidRapXiPlus(0),
    fHistPtGenCascMidRapXiPlus(0),
   
    fHistThetaGenCascXiPlus(0), 
    fHistPtFdblGenCascXiPlus(0),
    
    fHistThetaLambdaXiPlus(0), 
    fHistThetaBachXiPlus(0),
    
    fHistThetaMesDghterXiPlus(0), 
    fHistThetaBarDghterXiPlus(0),
    
    fHistPtBachXiPlus(0),
    fHistPtMesDghterXiPlus(0),
    fHistPtBarDghterXiPlus(0),
   
   // Omega-
   fHistEtaGenCascOmegaMinus(0),
   
    fHistYGenCascMidRapOmegaMinus(0),
    fHistEtaGenCascMidRapOmegaMinus(0),
    fHistThetaGenCascMidRapOmegaMinus(0),
    fHistPtGenCascMidRapOmegaMinus(0),
   
    fHistThetaGenCascOmegaMinus(0),
    fHistPtFdblGenCascOmegaMinus(0),
    
    fHistThetaLambdaOmegaMinus(0), 
    fHistThetaBachOmegaMinus(0),
    
    fHistThetaMesDghterOmegaMinus(0), 
    fHistThetaBarDghterOmegaMinus(0),
    
    fHistPtBachOmegaMinus(0),
    fHistPtMesDghterOmegaMinus(0),
    fHistPtBarDghterOmegaMinus(0),
   
   
   // Omega+
   fHistEtaGenCascOmegaPlus(0),
    
    fHistYGenCascMidRapOmegaPlus(0),
    fHistEtaGenCascMidRapOmegaPlus(0),
    fHistThetaGenCascMidRapOmegaPlus(0),
    fHistPtGenCascMidRapOmegaPlus(0),
   
    fHistThetaGenCascOmegaPlus(0),
    fHistPtFdblGenCascOmegaPlus(0),
    
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
	
	// - Effective mass histos for the cascade candidates associated with MC
    fHistAsMCMassXiMinus(0),		
    fHistAsMCMassXiPlus(0),		
    fHistAsMCMassOmegaMinus(0),
    fHistAsMCMassOmegaPlus(0),
	
	// - Generated Pt of the the cascade candidates associated with MC
    fHistAsMCGenPtXiMinus(0),		
    fHistAsMCGenPtXiPlus(0),		
    fHistAsMCGenPtOmegaMinus(0),
    fHistAsMCGenPtOmegaPlus(0),
	
	// - Generated Y of the the cascade candidates associated with MC
    fHistAsMCGenYXiMinus(0),		
    fHistAsMCGenYXiPlus(0),	
    fHistAsMCGenYOmegaMinus(0),
    fHistAsMCGenYOmegaPlus(0),
    
    	// - Generated Y Vs Generated Pt, for the cascade candidates associated with MC
    f2dHistAsMCGenYVsGenPtXiMinus(0),
    f2dHistAsMCGenYVsGenPtXiPlus(0),
    f2dHistAsMCGenYVsGenPtOmegaMinus(0),
    f2dHistAsMCGenYVsGenPtOmegaPlus(0),
    
    
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
    f2dHistAsMCResROmegaPlus(0)
	
    
		    
{
// Dummy constructor
}
     
       
     
     
//_____Non-default Constructor________________________________________________________________
AliAnalysisTaskCheckPerformanceCascade::AliAnalysisTaskCheckPerformanceCascade(const char *name) 
  : AliAnalysisTaskSE(name),
    fDebugCascade(0), fAnalysisType("ESD"), fCollidingSystems(0),
      
    	// - Cascade part initialisation
fListHistCascade(0),
    fHistMCTrackMultiplicity(0), 

    fHistEtaGenProton(0),
    fHistEtaGenAntiProton(0),
    
   // Xi-
    fHistEtaGenCascXiMinus(0),
      
    fHistYGenCascMidRapXiMinus(0),		
    fHistEtaGenCascMidRapXiMinus(0),
    fHistThetaGenCascMidRapXiMinus(0),
    fHistPtGenCascMidRapXiMinus(0),
   
    fHistThetaGenCascXiMinus(0),
    fHistPtFdblGenCascXiMinus(0),
    
    fHistThetaLambdaXiMinus(0), 
    fHistThetaBachXiMinus(0),
    
    fHistThetaMesDghterXiMinus(0), 
    fHistThetaBarDghterXiMinus(0),
    
    fHistPtBachXiMinus(0),
    fHistPtMesDghterXiMinus(0),
    fHistPtBarDghterXiMinus(0),
   
   
   // Xi+
    fHistEtaGenCascXiPlus(0),
      
    fHistYGenCascMidRapXiPlus(0),
    fHistEtaGenCascMidRapXiPlus(0),
    fHistThetaGenCascMidRapXiPlus(0),
    fHistPtGenCascMidRapXiPlus(0),
   
    fHistThetaGenCascXiPlus(0), 
    fHistPtFdblGenCascXiPlus(0),
    
    fHistThetaLambdaXiPlus(0), 
    fHistThetaBachXiPlus(0),
    
    fHistThetaMesDghterXiPlus(0), 
    fHistThetaBarDghterXiPlus(0),
    
    fHistPtBachXiPlus(0),
    fHistPtMesDghterXiPlus(0),
    fHistPtBarDghterXiPlus(0),
   
   // Omega-
    fHistEtaGenCascOmegaMinus(0),
   
    fHistYGenCascMidRapOmegaMinus(0),
    fHistEtaGenCascMidRapOmegaMinus(0),
    fHistThetaGenCascMidRapOmegaMinus(0),
    fHistPtGenCascMidRapOmegaMinus(0),
   
    fHistThetaGenCascOmegaMinus(0),
    fHistPtFdblGenCascOmegaMinus(0),
    
    fHistThetaLambdaOmegaMinus(0), 
    fHistThetaBachOmegaMinus(0),
    
    fHistThetaMesDghterOmegaMinus(0), 
    fHistThetaBarDghterOmegaMinus(0),
    
    fHistPtBachOmegaMinus(0),
    fHistPtMesDghterOmegaMinus(0),
    fHistPtBarDghterOmegaMinus(0),
   
   
   // Omega+
    fHistEtaGenCascOmegaPlus(0),
    
    fHistYGenCascMidRapOmegaPlus(0),
    fHistEtaGenCascMidRapOmegaPlus(0),
    fHistThetaGenCascMidRapOmegaPlus(0),
    fHistPtGenCascMidRapOmegaPlus(0),
   
    fHistThetaGenCascOmegaPlus(0),
    fHistPtFdblGenCascOmegaPlus(0),
    
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
	
	// - Effective mass histos for the cascade candidates associated with MC
    fHistAsMCMassXiMinus(0),		
    fHistAsMCMassXiPlus(0),		
    fHistAsMCMassOmegaMinus(0),
    fHistAsMCMassOmegaPlus(0),
	
	// - Generated Pt of the the cascade candidates associated with MC
    fHistAsMCGenPtXiMinus(0),		
    fHistAsMCGenPtXiPlus(0),		
    fHistAsMCGenPtOmegaMinus(0),
    fHistAsMCGenPtOmegaPlus(0),
	
	// - Generated Y of the the cascade candidates associated with MC
    fHistAsMCGenYXiMinus(0),		
    fHistAsMCGenYXiPlus(0),	
    fHistAsMCGenYOmegaMinus(0),
    fHistAsMCGenYOmegaPlus(0),
    
    	// - Generated Y Vs Generated Pt, for the cascade candidates associated with MC
    f2dHistAsMCGenYVsGenPtXiMinus(0),
    f2dHistAsMCGenYVsGenPtXiPlus(0),
    f2dHistAsMCGenYVsGenPtOmegaMinus(0),
    f2dHistAsMCGenYVsGenPtOmegaPlus(0),
    
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
    f2dHistAsMCResROmegaPlus(0)
	
    
    
  
    
   
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
  
  // - Plots for |y(MC)| < 1	
  if (!fHistYGenCascMidRapXiMinus) {
	  fHistYGenCascMidRapXiMinus = new TH1F("fHistYGenCascMidRapXiMinus", "Y of #Xi^{-} generated at mid rapidity;Y;Number of Casc", 240, -1.2, 1.2);
	  fListHistCascade->Add(fHistYGenCascMidRapXiMinus );
  }
  
  if (!fHistEtaGenCascMidRapXiMinus) {
	  fHistEtaGenCascMidRapXiMinus = new TH1F("fHistEtaGenCascMidRapXiMinus", "#eta of #Xi^{-} generated at mid rapidity;#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add(fHistEtaGenCascMidRapXiMinus );
  }
  
  if (!fHistThetaGenCascMidRapXiMinus) {
	  fHistThetaGenCascMidRapXiMinus = new TH1F("fHistThetaGenCascMidRapXiMinus", "#theta of #Xi^{-} generated at mid rapidity;#theta (deg);Number of Casc", 200, -10, 190);
	  fListHistCascade->Add(fHistThetaGenCascMidRapXiMinus );
  }
  
  if (!fHistPtGenCascMidRapXiMinus) {
	  fHistPtGenCascMidRapXiMinus = new TH1F("fHistPtGenCascMidRapXiMinus", "P_{t} of #Xi^{-} generated at mid rapidity;P_{t} (GeV/c);Number of Casc", 200, 0., 10.);
	  fListHistCascade->Add(fHistPtGenCascMidRapXiMinus );
  }
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascXiMinus) {
     fHistThetaGenCascXiMinus = new TH1F("fHistThetaGenCascXiMinus", "#theta of gen. #Xi^{-};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascXiMinus);
  }

  if (!fHistPtFdblGenCascXiMinus) {
	  fHistPtFdblGenCascXiMinus = new TH1F("fHistPtFdblGenCascXiMinus","P_{t} of findable generated #Xi^{-};P_{t} (GeV/c);Number of Casc", 200, 0, 10);
	  fListHistCascade->Add(fHistPtFdblGenCascXiMinus );
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
  
  // - Plots for |y(MC)| < 1	
  if (!fHistYGenCascMidRapXiPlus) {
	  fHistYGenCascMidRapXiPlus = new TH1F("fHistYGenCascMidRapXiPlus", "Y of #Xi^{+} generated at mid rapidity;Y;Number of Casc", 240, -1.2, 1.2);
	  fListHistCascade->Add(fHistYGenCascMidRapXiPlus );
  }
  
  if (!fHistEtaGenCascMidRapXiPlus) {
	  fHistEtaGenCascMidRapXiPlus = new TH1F("fHistEtaGenCascMidRapXiPlus", "#eta of #Xi^{+} generated at mid rapidity;#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add(fHistEtaGenCascMidRapXiPlus );
  }
  
  if (!fHistThetaGenCascMidRapXiPlus) {
	  fHistThetaGenCascMidRapXiPlus = new TH1F("fHistThetaGenCascMidRapXiPlus", "#theta of #Xi^{+} generated at mid rapidity;#theta (deg);Number of Casc", 200, -10, 190);
	  fListHistCascade->Add(fHistThetaGenCascMidRapXiPlus );
  }
  
  if (!fHistPtGenCascMidRapXiPlus) {
	  fHistPtGenCascMidRapXiPlus = new TH1F("fHistPtGenCascMidRapXiPlus", "P_{t} of #Xi^{+} generated at mid rapidity;P_{t} (GeV/c);Number of Casc", 200, 0., 10.);
	  fListHistCascade->Add(fHistPtGenCascMidRapXiPlus );
  }
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascXiPlus) {
     fHistThetaGenCascXiPlus = new TH1F("fHistThetaGenCascXiPlus", "#theta of gen. #Xi^{+};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascXiPlus);
  }
 
  if (!fHistPtFdblGenCascXiPlus) {
	  fHistPtFdblGenCascXiPlus = new TH1F("fHistPtFdblGenCascXiPlus","P_{t} of findable generated #Xi^{+};P_{t} (GeV/c);Number of Casc", 200, 0, 10);
	  fListHistCascade->Add(fHistPtFdblGenCascXiPlus );
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
  
  // - Plots for |y(MC)| < 1	
  if (!fHistYGenCascMidRapOmegaMinus) {
	  fHistYGenCascMidRapOmegaMinus = new TH1F("fHistYGenCascMidRapOmegaMinus", "Y of #Omega^{-} generated at mid rapidity;Y;Number of Casc", 240, -1.2, 1.2);
	  fListHistCascade->Add(fHistYGenCascMidRapOmegaMinus );
  }
  
  if (!fHistEtaGenCascMidRapOmegaMinus) {
	  fHistEtaGenCascMidRapOmegaMinus = new TH1F("fHistEtaGenCascMidRapOmegaMinus", "#eta of #Omega^{-} generated at mid rapidity;#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add(fHistEtaGenCascMidRapOmegaMinus );
  }
  
  if (!fHistThetaGenCascMidRapOmegaMinus) {
	  fHistThetaGenCascMidRapOmegaMinus = new TH1F("fHistThetaGenCascMidRapOmegaMinus", "#theta of #Omega^{-} generated at mid rapidity;#theta (deg);Number of Casc", 200, -10, 190);
	  fListHistCascade->Add(fHistThetaGenCascMidRapOmegaMinus );
  }
  
  if (!fHistPtGenCascMidRapOmegaMinus) {
	  fHistPtGenCascMidRapOmegaMinus = new TH1F("fHistPtGenCascMidRapOmegaMinus", "P_{t} of #Omega^{-} generated at mid rapidity;P_{t} (GeV/c);Number of Casc", 200, 0., 10.);
	  fListHistCascade->Add(fHistPtGenCascMidRapOmegaMinus );
  }
  
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascOmegaMinus) {
     fHistThetaGenCascOmegaMinus = new TH1F("fHistThetaGenCascOmegaMinus", "#theta of gen. #Omega^{-};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascOmegaMinus);
  }
 
  if (!fHistPtFdblGenCascOmegaMinus) {
	  fHistPtFdblGenCascOmegaMinus = new TH1F("fHistPtFdblGenCascOmegaMinus","P_{t} of findable generated #Omega^{-};P_{t} (GeV/c);Number of Casc", 200, 0, 10);
	  fListHistCascade->Add(fHistPtFdblGenCascOmegaMinus );
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
  
  // - Plots for |y(MC)| < 1	
  if (!fHistYGenCascMidRapOmegaPlus) {
	  fHistYGenCascMidRapOmegaPlus = new TH1F("fHistYGenCascMidRapOmegaPlus", "Y of #Omega^{+} generated at mid rapidity;Y;Number of Casc", 240, -1.2, 1.2);
	  fListHistCascade->Add(fHistYGenCascMidRapOmegaPlus );
  }
  
  if (!fHistEtaGenCascMidRapOmegaPlus) {
	  fHistEtaGenCascMidRapOmegaPlus = new TH1F("fHistEtaGenCascMidRapOmegaPlus", "#eta of #Omega^{+} generated at mid rapidity;#eta;Number of Casc", 100, -5, 5);
	  fListHistCascade->Add(fHistEtaGenCascMidRapOmegaPlus );
  }
  
  if (!fHistThetaGenCascMidRapOmegaPlus) {
	  fHistThetaGenCascMidRapOmegaPlus = new TH1F("fHistThetaGenCascMidRapOmegaPlus", "#theta of #Omega^{+} generated at mid rapidity;#theta (deg);Number of Casc", 200, -10, 190);
	  fListHistCascade->Add(fHistThetaGenCascMidRapOmegaPlus );
  }
  
  if (!fHistPtGenCascMidRapOmegaPlus) {
	  fHistPtGenCascMidRapOmegaPlus = new TH1F("fHistPtGenCascMidRapOmegaPlus", "P_{t} of #Omega^{+} generated at mid rapidity;P_{t} (GeV/c);Number of Casc", 200, 0., 10.);
	  fListHistCascade->Add(fHistPtGenCascMidRapOmegaPlus );
  }
  
  
  
  // - Info at the generation level of multi-strange particle
  
  if (!fHistThetaGenCascOmegaPlus) {
     fHistThetaGenCascOmegaPlus = new TH1F("fHistThetaGenCascOmegaPlus", "#theta of gen. #Omega^{+};#theta;Number of Casc.", 200, -10, 190);
     fListHistCascade->Add(fHistThetaGenCascOmegaPlus);
  }
 
  if (!fHistPtFdblGenCascOmegaPlus) {
	  fHistPtFdblGenCascOmegaPlus = new TH1F("fHistPtFdblGenCascOmegaPlus","P_{t} of findable generated #Omega^{+};P_{t} (GeV/c);Number of Casc", 200, 0, 10);
	  fListHistCascade->Add(fHistPtFdblGenCascOmegaPlus );
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
  
  
  		// - Generated Pt of the the cascade candidates associated with MC
  
  if (!fHistAsMCGenPtXiMinus) {
	  fHistAsMCGenPtXiMinus = new TH1F("fHistAsMCGenPtXiMinus", "P_{t} of gen. #Xi^{-} (associated);P_{t} (GeV/c);Number of Casc", 200, 0., 10.);
	  fListHistCascade->Add(fHistAsMCGenPtXiMinus );
  }
  
  if (!fHistAsMCGenPtXiPlus) {
	  fHistAsMCGenPtXiPlus = new TH1F("fHistAsMCGenPtXiPlus", "P_{t} of gen. #Xi^{+} (associated);P_{t} (GeV/c);Number of Casc", 200, 0., 10.);
	  fListHistCascade->Add(fHistAsMCGenPtXiPlus );
  }
  
  if (!fHistAsMCGenPtOmegaMinus) {
	  fHistAsMCGenPtOmegaMinus = new TH1F("fHistAsMCGenPtOmegaMinus", "P_{t} of gen. #Omega^{-} (associated);P_{t} (GeV/c);Number of Casc", 200, 0., 10.);
	  fListHistCascade->Add(fHistAsMCGenPtOmegaMinus );
  }
  
  if (!fHistAsMCGenPtOmegaPlus) {
	  fHistAsMCGenPtOmegaPlus = new TH1F("fHistAsMCGenPtOmegaPlus", "P_{t} of gen. #Omega^{+} (associated);P_{t} (GeV/c);Number of Casc", 200, 0., 10.);
	  fListHistCascade->Add(fHistAsMCGenPtOmegaPlus );
  }
  
  
  		// - Generated Y of the the cascade candidates associated with MC
  
  if (!fHistAsMCGenYXiMinus) {
	  fHistAsMCGenYXiMinus = new TH1F("fHistAsMCGenYXiMinus", "Rapidity, Y of gen. #Xi^{-} (associated);Rapidity, Y;Number of Casc", 240, -1.2, 1.2);
	  fListHistCascade->Add(fHistAsMCGenYXiMinus );
  }
  
  if (!fHistAsMCGenYXiPlus) {
	  fHistAsMCGenYXiPlus = new TH1F("fHistAsMCGenYXiPlus", "Rapidity, Y of gen. #Xi^{+} (associated);Rapidity, Y;Number of Casc", 240, -1.2, 1.2);
	  fListHistCascade->Add(fHistAsMCGenYXiPlus );
  }
  
  if (!fHistAsMCGenYOmegaMinus) {
	  fHistAsMCGenYOmegaMinus = new TH1F("fHistAsMCGenYOmegaMinus", "Rapidity, Y of gen. #Omega^{-} (associated);Rapidity, Y;Number of Casc", 240, -1.2, 1.2);
	  fListHistCascade->Add(fHistAsMCGenYOmegaMinus );
  }
  
  if (!fHistAsMCGenYOmegaPlus) {
	  fHistAsMCGenYOmegaPlus = new TH1F("fHistAsMCGenYOmegaPlus", "Rapidity, Y of gen. #Omega^{+} (associated);Rapidity, Y;Number of Casc", 240, -1.2, 1.2);
	  fListHistCascade->Add(fHistAsMCGenYOmegaPlus );
  }
  
  		// - Generated Y Vs Generated Pt, for the cascade candidates associated with MC
  
  if (!f2dHistAsMCGenYVsGenPtXiMinus) {
  	f2dHistAsMCGenYVsGenPtXiMinus = new TH2F("f2dHistAsMCGenYVsGenPtXiMinus", "Y Vs P_{t} of gen. #Xi^{-} (associated);Rapidity, Y;P_{t} (GeV/c)",240, -1.2, 1.2, 200, 0., 10.);
	  fListHistCascade->Add(f2dHistAsMCGenYVsGenPtXiMinus );
  }
  
  if (!f2dHistAsMCGenYVsGenPtXiPlus) {
	  f2dHistAsMCGenYVsGenPtXiPlus = new TH2F("f2dHistAsMCGenYVsGenPtXiPlus", "Y Vs P_{t} of gen. #Xi^{+} (associated);Rapidity, Y;P_{t} (GeV/c)",240, -1.2, 1.2, 200, 0., 10.);
	  fListHistCascade->Add(f2dHistAsMCGenYVsGenPtXiPlus );
  }
  
  if (!f2dHistAsMCGenYVsGenPtOmegaMinus) {
	  f2dHistAsMCGenYVsGenPtOmegaMinus = new TH2F("f2dHistAsMCGenYVsGenPtOmegaMinus", "Y Vs P_{t} of gen. #Omega^{-} (associated);Rapidity, Y;P_{t} (GeV/c)",240, -1.2, 1.2, 200, 0., 10.);
	  fListHistCascade->Add(f2dHistAsMCGenYVsGenPtOmegaMinus );
  }
  
  if (!f2dHistAsMCGenYVsGenPtOmegaPlus) {
	  f2dHistAsMCGenYVsGenPtOmegaPlus = new TH2F("f2dHistAsMCGenYVsGenPtOmegaPlus", "Y Vs P_{t} of gen. #Omega^{+} (associated);Rapidity, Y;P_{t} (GeV/c)",240, -1.2, 1.2, 200, 0., 10.);
	  fListHistCascade->Add(f2dHistAsMCGenYVsGenPtOmegaPlus );
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
	  f2dHistAsMCResPtXiMinus = new TH2F( "f2dHistAsMCResPtXiMinus", "Resolution in Pt reconstruction for #Xi^{-}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC} (%)", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtXiMinus);
  }
  
  if(! f2dHistAsMCResPtXiPlus) {
	  f2dHistAsMCResPtXiPlus = new TH2F( "f2dHistAsMCResPtXiPlus", "Resolution in Pt reconstruction for #Xi^{+}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC} (%)", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtXiPlus);
  }
  
  if(! f2dHistAsMCResPtOmegaMinus) {
	  f2dHistAsMCResPtOmegaMinus = new TH2F( "f2dHistAsMCResPtOmegaMinus", "Resolution in Pt reconstruction for #Omega^{-}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC} (%)", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtOmegaMinus);
  }
  
  if(! f2dHistAsMCResPtOmegaPlus) {
	  f2dHistAsMCResPtOmegaPlus = new TH2F( "f2dHistAsMCResPtOmegaPlus", "Resolution in Pt reconstruction for #Omega^{+}; Pt_{MC} (GeV/c); (Pt_{reco} - Pt_{MC}) / Pt_{MC} (%)", 200, 0., 10., 200, -0.1, 0.1);
	  fListHistCascade->Add(f2dHistAsMCResPtOmegaPlus);
  }
  
  		// - Resolution in R(2D) as function of generated R
  
  if(! f2dHistAsMCResRXiMinus) {
	  f2dHistAsMCResRXiMinus = new TH2F( "f2dHistAsMCResRXiMinus", "Resolution in transv. position for #Xi^{-}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC} (%)", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResRXiMinus);
  }
  
  if(! f2dHistAsMCResRXiPlus) {
	  f2dHistAsMCResRXiPlus = new TH2F( "f2dHistAsMCResRXiPlus", "Resolution in transv. position for #Xi^{+}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC} (%)", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResRXiPlus);
  }
  
  if(! f2dHistAsMCResROmegaMinus) {
	  f2dHistAsMCResROmegaMinus = new TH2F( "f2dHistAsMCResROmegaMinus", "Resolution in transv. position for #Omega^{-}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC} (%)", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResROmegaMinus);
  }
  
  if(! f2dHistAsMCResROmegaPlus) {
	  f2dHistAsMCResROmegaPlus = new TH2F( "f2dHistAsMCResROmegaPlus", "Resolution in transv. position for #Omega^{+}; R_{MC} (cm); (R_{reco} - R_{MC}) / R_{MC} (%)", 450, 0., 45.0, 240, -0.3, 0.3);
	  fListHistCascade->Add(f2dHistAsMCResROmegaPlus);
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
	

  //	cout << "Name of the accessed file :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;

  //	cout << "Tree characteristics ..." << endl;
  //	fInputHandler->GetTree()->Print("toponly");
  //	fInputHandler->GetTree()->GetBranch("PrimaryVertex")->Print();
  //	fInputHandler->GetTree()->GetBranch("SPDVertex")->Print();

 
 
 
 

 
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
	
	TH1F *lHistYGenCascMidRap     = 0;
	TH1F *lHistEtaGenCascMidRap   = 0;
	TH1F *lHistThetaGenCascMidRap = 0;
	TH1F *lHistPtGenCascMidRap    = 0;
	
	TH1F *lHistThetaGenCasc       = 0;
	TH1F *lHistPtFdblGenCasc      = 0;
	TH1F *lHistThetaLambda        = 0;
	TH1F *lHistThetaBach          = 0;
	TH1F *lHistThetaBarDghter     = 0;
	TH1F *lHistThetaMesDghter     = 0;
	TH1F *lHistPtBach             = 0;
	TH1F *lHistPtBarDghter        = 0;
	TH1F *lHistPtMesDghter        = 0;


for(Int_t CascType = 1; CascType < 5; CascType++)
{
       
switch (CascType)
  {
    case 1: // Xi-
         lPdgCodeCasc       =   3312;  //Xi-
         lPdgCodeBach       =   -211;  //Pi-
         lPdgCodeLambda     =   3122;  //Lambda0
         lPdgCodeDghtMesV0  =   -211;  //Pi-
         lPdgCodeDghtBarV0  =   2212;  //Proton 
	 	
	 	// any Xi-
	 lHistEtaGenCasc     = (TH1F*)fHistEtaGenCascXiMinus;
	 
	 	// cascades within |y| < 1
	 lHistYGenCascMidRap     	= fHistYGenCascMidRapXiMinus;		
	 lHistEtaGenCascMidRap		= fHistEtaGenCascMidRapXiMinus;
	 lHistThetaGenCascMidRap	= fHistThetaGenCascMidRapXiMinus;	
	 lHistPtGenCascMidRap		= fHistPtGenCascMidRapXiMinus;		
	 	
	 	// cascades generated within acceptance (cut in pt + theta)
	 lHistThetaGenCasc   = fHistThetaGenCascXiMinus;
	 lHistPtFdblGenCasc  = fHistPtFdblGenCascXiMinus;
	 lHistThetaLambda    = fHistThetaLambdaXiMinus;
	 lHistThetaBach      = fHistThetaBachXiMinus;
	 lHistThetaBarDghter = fHistThetaBarDghterXiMinus;
	 lHistThetaMesDghter = fHistThetaMesDghterXiMinus;
	 lHistPtBach	     = fHistPtBachXiMinus;
	 lHistPtBarDghter    = fHistPtBarDghterXiMinus;
	 lHistPtMesDghter    = fHistPtMesDghterXiMinus;
        break; 
           
    case 2: // Xi+
         lPdgCodeCasc        =  -3312;  //Xi+
         lPdgCodeBach        =    211;  //Pi+
         lPdgCodeLambda      =  -3122;  //AntiLambda0
         lPdgCodeDghtMesV0   =    211;  //Pi+
         lPdgCodeDghtBarV0   =  -2212;  //AntiProton  
	 
	 	// any Xi+
	 lHistEtaGenCasc     		= fHistEtaGenCascXiPlus;
	 
	 	// cascades within |y| < 1
	 lHistYGenCascMidRap     	= fHistYGenCascMidRapXiPlus;		
	 lHistEtaGenCascMidRap		= fHistEtaGenCascMidRapXiPlus;
	 lHistThetaGenCascMidRap	= fHistThetaGenCascMidRapXiPlus;	
	 lHistPtGenCascMidRap		= fHistPtGenCascMidRapXiPlus;		
	 	
	 	// cascades generated within acceptance (cut in pt + theta)	 
	 lHistThetaGenCasc   = fHistThetaGenCascXiPlus;
	 lHistPtFdblGenCasc  = fHistPtFdblGenCascXiPlus;
	 lHistThetaLambda    = fHistThetaLambdaXiPlus;
	 lHistThetaBach      = fHistThetaBachXiPlus;
	 lHistThetaBarDghter = fHistThetaBarDghterXiPlus;
	 lHistThetaMesDghter = fHistThetaMesDghterXiPlus;
	 lHistPtBach	     = fHistPtBachXiPlus;
	 lHistPtBarDghter    = fHistPtBarDghterXiPlus;
	 lHistPtMesDghter    = fHistPtMesDghterXiPlus;  
    	break;
   
    case 3: // Omega-
    	 lPdgCodeCasc       =   3334;  //Omega-
         lPdgCodeBach       =   -321;  //K-
         lPdgCodeLambda     =   3122;  //Lambda0
         lPdgCodeDghtMesV0  =   -211;  //Pi-
         lPdgCodeDghtBarV0  =   2212;  //Proton 
	 
	 	// any Omega-
	 lHistEtaGenCasc     = fHistEtaGenCascOmegaMinus;
	 	
	 	// cascades within |y| < 1
	 lHistYGenCascMidRap     	= fHistYGenCascMidRapOmegaMinus;		
	 lHistEtaGenCascMidRap		= fHistEtaGenCascMidRapOmegaMinus;
	 lHistThetaGenCascMidRap	= fHistThetaGenCascMidRapOmegaMinus;	
	 lHistPtGenCascMidRap		= fHistPtGenCascMidRapOmegaMinus;		
	 	
	 	// cascades generated within acceptance (cut in pt + theta)
	 lHistThetaGenCasc   = fHistThetaGenCascOmegaMinus;
	 lHistPtFdblGenCasc  = fHistPtFdblGenCascOmegaMinus;
	 lHistThetaLambda    = fHistThetaLambdaOmegaMinus;
	 lHistThetaBach      = fHistThetaBachOmegaMinus;
	 lHistThetaBarDghter = fHistThetaBarDghterOmegaMinus;
	 lHistThetaMesDghter = fHistThetaMesDghterOmegaMinus;
	 lHistPtBach	     = fHistPtBachOmegaMinus;
	 lHistPtBarDghter    = fHistPtBarDghterOmegaMinus;
	 lHistPtMesDghter    = fHistPtMesDghterOmegaMinus;   
        break;
    
    case 4:  // Omega+
         lPdgCodeCasc       =  -3334;  //Omega+
         lPdgCodeBach       =    321;  //K+
         lPdgCodeLambda     =  -3122;  //AntiLambda0
         lPdgCodeDghtMesV0  =    211;  //Pi+
         lPdgCodeDghtBarV0  =  -2212;  //AntiProton 
	 
	 	// any Omega+
	 lHistEtaGenCasc     = fHistEtaGenCascOmegaPlus;
	 
	 	// cascades within |y| < 1
	 lHistYGenCascMidRap     	= fHistYGenCascMidRapOmegaPlus;		
	 lHistEtaGenCascMidRap		= fHistEtaGenCascMidRapOmegaPlus;
	 lHistThetaGenCascMidRap	= fHistThetaGenCascMidRapOmegaPlus;	
	 lHistPtGenCascMidRap		= fHistPtGenCascMidRapOmegaPlus;		
	 	
	 	// cascades generated within acceptance (cut in pt + theta)
	 lHistThetaGenCasc   = fHistThetaGenCascOmegaPlus;
	 lHistPtFdblGenCasc  = fHistPtFdblGenCascOmegaPlus;
	 lHistThetaLambda    = fHistThetaLambdaOmegaPlus;
	 lHistThetaBach      = fHistThetaBachOmegaPlus;
	 lHistThetaBarDghter = fHistThetaBarDghterOmegaPlus;
	 lHistThetaMesDghter = fHistThetaMesDghterOmegaPlus;
	 lHistPtBach	     = fHistPtBachOmegaPlus;
	 lHistPtBarDghter    = fHistPtBarDghterOmegaPlus;
	 lHistPtMesDghter    = fHistPtMesDghterOmegaPlus;  
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
   		//cout << "Xi- ds boucle " << iCurrentLabelStack << "/ " << iNumberOfPrimaries << endl;
		
		// -  Xi level ... _____________________________________________________________
		TParticle* xiMC = 0x0;
			   xiMC = lCurrentParticle;
		if(!xiMC){
			Printf("MC TParticle pointer to Cascade = 0x0 ! Skip ...");
			continue;
		
		}
		
		Double_t lRadToDeg = 180.0/TMath::Pi();
		
		// Fill the histo in pseudo-rapidity : = any generated Xi, not necessarily within the acceptance
		lHistEtaGenCasc->Fill( xiMC->Eta() );
		
	 
	 	// cascades within |y| < 1
		Double_t lRapXi = 0.5*TMath::Log((xiMC->Energy() + xiMC->Pz()) / (xiMC->Energy() - xiMC->Pz() +1.e-13));
		
		if( TMath::Abs(lRapXi) < 1) {
			lHistYGenCascMidRap 	->Fill( lRapXi    			);    	
			lHistEtaGenCascMidRap	->Fill( xiMC->Eta() 			); 
			lHistThetaGenCascMidRap ->Fill( lRadToDeg * xiMC->Theta()	);
			lHistPtGenCascMidRap	->Fill( xiMC->Pt() 			);
		}
		
		// Check the emission of particle stays within the acceptance of the detector (cut in pt + theta)
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
		
			// Check the emission of particle stays within the acceptance of the detector
			if( lLambda->Theta() < TMath::Pi()/4.0  ||    lLambda->Theta() > 3.0*TMath::Pi()/4.0 ) continue;
			if( lBach->Theta() < TMath::Pi()/4.0    ||    lBach->Theta() > 3.0*TMath::Pi()/4.0 )   continue;
		
			if( lBach->Pt() < 0.2 ) continue; 
			
		
		
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
			
			if( lDghtBarV0->Pt() < 0.5 ) continue;
			if( lDghtMesV0->Pt() < 0.2 ) continue;
			
			
			
		// - Just to know which file is currently open : locate the file containing Xi 
		//cout << "Name of the file containing generated Xi :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() 
		//						     <<  endl;	
			
		// - Filling histos ... _________________________________________________________________	
			lHistThetaGenCasc	->Fill( lRadToDeg * xiMC->Theta()  );
			lHistPtFdblGenCasc	->Fill( xiMC->Pt() );
			
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
			
			// - Re-initialisation of the local TH1F pointers
			lHistEtaGenCasc         = 0x0;
			
			lHistYGenCascMidRap     = 0x0;
			lHistEtaGenCascMidRap   = 0x0;
			lHistThetaGenCascMidRap = 0x0;
			lHistPtGenCascMidRap	= 0x0;
			
			lHistThetaGenCasc       = 0x0;
			lHistPtFdblGenCasc      = 0x0;
			lHistThetaLambda        = 0x0;
			lHistThetaBach          = 0x0;
			lHistThetaBarDghter     = 0x0;
			lHistThetaMesDghter     = 0x0;
			lHistPtBach	        = 0x0;
			lHistPtBarDghter        = 0x0;
			lHistPtMesDghter        = 0x0;	
			
			
		}// end if current particle = Xi-
	     
     }// This is the end of the loop on primaries

} // end of loop over the different types of cascades (Xi-+, Omega-+)
 	
 
 
//__________________________________________________________________________	
// Part 2 - Loop over the reconstructed candidates
  
  
// Temporary way : AOD awareness of the code to be developed  
if(fAnalysisType == "AOD") return;


for (Int_t iXi = 0; iXi < ncascades; iXi++) 
{// This is the begining of the Cascade loop
		
	AliESDcascade *xiESD = lESDevent->GetCascade(iXi);
          
	if (!xiESD) continue;
	
	// -  Step 1 : Preparing the general info about of the event
	//-------------
	//const AliESDVertex *lPrimaryVtx = lESDevent->GetPrimaryVertex();  // get the best vtx available
	//	Double_t lPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
	//		 lPrimaryVtx->GetXYZ( lPrimaryVtxPos );
	
	// Double_t lMagneticField = lESDevent->GetMagneticField( );
	
	
	// - Step 2 : Connection to daughter tracks of the current cascade
	//-------------
			
		UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xiESD->GetPindex() );
		UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xiESD->GetNindex() );
		UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xiESD->GetBindex() );
		// abs value not needed ; the index should always be positive (!= label ...)
      
	AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
	AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
	AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx  );
	if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
		Printf("ERROR: Could not retrieve one of the 3 daughter tracks of the cascade ...");
		continue;
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
	
	
	
	
	// - Step 4 : MC association
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
	
		if(lblMotherPosV0Dghter != lblMotherNegV0Dghter) continue; // same mother
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

	
	// - Step 5 : Plots around the cascade candidates associated with MC
	//-------------	
	
	Double_t lmcPt             = mcMotherBach->Pt();
	Double_t lmcRapXi          = 0.5*TMath::Log( (mcMotherBach->Energy() + mcMotherBach->Pz()) / 
						     (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) );
	Double_t lmcEta            = mcMotherBach->Eta();
	Double_t lmcTransvRadius   = mcBach->R(); // to get the decay point of Xi, = the production vertex of Bachelor ...
	
	Double_t lrecoPt           = xiESD->Pt();
	Double_t lrecoTransvRadius = TMath::Sqrt( xiESD->Xv() * xiESD->Xv() + xiESD->Yv() * xiESD->Yv() );
	
	// - Effective mass histos for the cascade candidates associated with MC
	
	if( lChargeXi < 0 && lAssoXiMinus){	
		fHistAsMCMassXiMinus	      ->Fill( lInvMassXiMinus  );
		fHistAsMCGenPtXiMinus         ->Fill( lmcPt            );
		fHistAsMCGenYXiMinus          ->Fill( lmcRapXi         );
		f2dHistAsMCGenYVsGenPtXiMinus ->Fill( lmcRapXi, lmcPt  );
		fHistAsMCGenEtaXiMinus        ->Fill( lmcEta           );
		f2dHistAsMCResPtXiMinus       ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResRXiMinus        ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
	}
	
	else if( lChargeXi > 0 && lAssoXiPlus){	
		fHistAsMCMassXiPlus	      ->Fill( lInvMassXiPlus   );
		fHistAsMCGenPtXiPlus          ->Fill( lmcPt            );
		fHistAsMCGenYXiPlus           ->Fill( lmcRapXi         );
		f2dHistAsMCGenYVsGenPtXiPlus  ->Fill( lmcRapXi, lmcPt  );
		fHistAsMCGenEtaXiPlus         ->Fill( lmcEta           );
		f2dHistAsMCResPtXiPlus        ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResRXiPlus         ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
	}
	
	else if( lChargeXi < 0 && lAssoOmegaMinus){	
		fHistAsMCMassOmegaMinus          ->Fill( lInvMassOmegaMinus );
		fHistAsMCGenPtOmegaMinus         ->Fill( lmcPt              );
		fHistAsMCGenYOmegaMinus          ->Fill( lmcRapXi           );
		f2dHistAsMCGenYVsGenPtOmegaMinus ->Fill( lmcRapXi, lmcPt    );
		fHistAsMCGenEtaOmegaMinus        ->Fill( lmcEta             );
		f2dHistAsMCResPtOmegaMinus       ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResROmegaMinus        ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
	}
	
	else if( lChargeXi > 0 && lAssoOmegaPlus){	
		fHistAsMCMassOmegaPlus           ->Fill( lInvMassOmegaPlus );
		fHistAsMCGenPtOmegaPlus          ->Fill( lmcPt             );
		fHistAsMCGenYOmegaPlus           ->Fill( lmcRapXi          );
		f2dHistAsMCGenYVsGenPtOmegaPlus  ->Fill( lmcRapXi, lmcPt   );
		fHistAsMCGenEtaOmegaPlus         ->Fill( lmcEta            );
		f2dHistAsMCResPtOmegaPlus        ->Fill( lmcPt,           (lrecoPt - lmcPt)/ lmcPt );
		f2dHistAsMCResROmegaPlus         ->Fill( lmcTransvRadius, (lrecoTransvRadius - lmcTransvRadius)/ lmcTransvRadius    );
	}
	
	
}// End of loop over reconstructed cascades
 
 
 
 
  // Post output data.
 PostData(1, fListHistCascade);
}      










//________________________________________________________________________
void AliAnalysisTaskCheckPerformanceCascade::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fHistMCTrackMultiplicity = dynamic_cast<TH1F*> (  ((TList*)GetOutputData(1))->FindObject("fHistMCTrackMultiplicity")  );
  if (!fHistMCTrackMultiplicity) {
    Printf("ERROR: fHistMCTrackMultiplicity not available");
    return;
  }
  
   
  TCanvas *c2 = new TCanvas("AliAnalysisTaskCheckPerformanceCascade","Multiplicity",10,10,510,510);
  c2->cd(1)->SetLogy();

  fHistMCTrackMultiplicity->SetMarkerStyle(22);
  fHistMCTrackMultiplicity->DrawCopy("E");
 // fHistV0Multiplicity->SetMarkerStyle(26);
 // fHistV0Multiplicity->DrawCopy("ESAME");

}
