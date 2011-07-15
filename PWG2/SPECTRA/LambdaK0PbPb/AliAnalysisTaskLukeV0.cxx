 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

/* $Id: AliAnalysisTaskLukeV0.cxx 46301 2011-01-06 14:25:27Z agheata $ */

/* AliAnalysisTaskLukeV0.cxx
 *
 * Task analysing lambda, antilambda & K0 spectra
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 * Adapted by Luke Hanratty
 *
 */

#include "AliAnalysisTaskLukeV0.h"

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TPDGCode.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliPID.h"
#include "AliPIDResponse.h"

ClassImp(AliAnalysisTaskLukeV0)

//________________________________________________________________________
AliAnalysisTaskLukeV0::AliAnalysisTaskLukeV0() // All data members should be initialised here
   :AliAnalysisTaskSE(),
    fOutputList(0),
    fTrackCuts(0),
	fPIDResponse(0),
	fHistCosPA(0), 
	fHistDCAV0Daughters(0), 
	fHistDecayL(0), 
	fHistImpactxyN(0), 
	fHistImpactzN(0), 
	fHistImpactxyP(0), 
	fHistImpactzP(0), 
	fHistMCLambdacTau(0),
	fHistMCLambdaNotcTau(0),
	fHistMCLambdaDecayL(0),
	fHistMCLambdaNotDecayL(0),
	fHistMcNLambdaPrimary(0),
	fHistMcNLambda(0),
	fHistMcNAntilambda(0),
	fHistMcNKshort(0),
	fHistMK0(0), 
	fHistMLa(0), 
	fHistMLb(0), 
	fHistNLambda(0), 
	fHistNV0(0), 
	fHistPtV0(0), 
	fHistPVZ(0), 
	fHistTauLa(0), 
	fHistV0Z(0), 
	fHistZ(0),
	fHistBetheBlochTPCNeg(0), 
	fHistBetheBlochTPCPos(0), 
	fHistImpactxyImpactz(0), 
	fHistMcPMK0Pt(0),
	fHistMcPMLaPt(0),
	fHistMcPMLbPt(0),
	fHistMcV0MK0Pt(0),
	fHistMcV0MLaPt(0),
	fHistMcV0MLbPt(0),
	fHistMK0Pt(0), 
	fHistMLaPt(0), 
	fHistMLbPt(0), 
	fHistPtArm(0), 
	fHistRZ(0), 
	fHistXZ(0), 
	fHistYZ(0) // The last in the above list should not have a comma after it
{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskLukeV0::AliAnalysisTaskLukeV0(const char *name) // All data members should be initialised here
   :AliAnalysisTaskSE(name),
	fOutputList(0),
	fTrackCuts(0),
	fPIDResponse(0),
	fHistCosPA(0), 
	fHistDCAV0Daughters(0), 
	fHistDecayL(0), 
	fHistImpactxyN(0), 
	fHistImpactzN(0), 
	fHistImpactxyP(0), 
	fHistImpactzP(0), 
	fHistMCLambdacTau(0),
	fHistMCLambdaNotcTau(0),
	fHistMCLambdaDecayL(0),
	fHistMCLambdaNotDecayL(0),
	fHistMcNLambdaPrimary(0),
	fHistMcNLambda(0),
	fHistMcNAntilambda(0),
	fHistMcNKshort(0),
	fHistMK0(0), 
	fHistMLa(0), 
	fHistMLb(0), 
	fHistNLambda(0), 
	fHistNV0(0), 
	fHistPtV0(0), 
	fHistPVZ(0), 
	fHistTauLa(0), 
	fHistV0Z(0), 
	fHistZ(0),
	fHistBetheBlochTPCNeg(0), 
	fHistBetheBlochTPCPos(0), 
	fHistImpactxyImpactz(0), 
	fHistMcPMK0Pt(0),
	fHistMcPMLaPt(0),
	fHistMcPMLbPt(0),
	fHistMcV0MK0Pt(0),
	fHistMcV0MLaPt(0),
	fHistMcV0MLbPt(0),
	fHistMK0Pt(0), 
	fHistMLaPt(0), 
	fHistMLbPt(0), 
	fHistPtArm(0), 
	fHistRZ(0), 
	fHistXZ(0), 
	fHistYZ(0) // The last in the above list should not have a comma after it
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskLukeV0::~AliAnalysisTaskLukeV0()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fOutputList;
    }
    if (fTrackCuts) delete fTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskLukeV0::UserCreateOutputObjects()
{
    // Create histograms
    // Called once (on the worker node)
        
    fOutputList = new TList();
    fOutputList->SetOwner();  // IMPORTANT!
    
	fTrackCuts = new AliESDtrackCuts();	
	
	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
	fPIDResponse = inputHandler->GetPIDResponse();

	// lambda plot parameters
	int div = 96;
	float max = 1.2;
	float min = 1.08;
	
	// Create remaining histograms
	// TH1F first
	fHistCosPA = new	TH1F("fHistCosPA", "Cosine of Pointing Angle of V0s; Cos PA; N(v0s)",202,0.8,1.01);
	fHistDCAV0Daughters = new	TH1F("fHistDCAV0Daughters", "DCA between V0 daughters; DCA (cm); N V0s", 100, 0, 2);
	fHistDecayL = new	TH1F("fHistDecayL", "Distance between V0 and PV; Distance(cm); N(v0s)",200,-0.1,30);
	fHistImpactxyN = new	TH1F("fHistImpactxyN", "RSM DCA between negative particle and primary vertex in xy plane; RSM DCA (cm); N(v0s)",100,0,1);
	fHistImpactzN = new	TH1F("fHistImpactzN", "RSM DCA between negative particle and primary vertex in z direction; RSM DCA (cm); N(v0s)",100,0,1);
	fHistImpactxyP = new	TH1F("fHistImpactxyP", "RSM DCA between positive particle and primary vertex in xy plane; RSM DCA (cm); N(v0s)",100,0,1);
	fHistImpactzP = new	TH1F("fHistImpactzP", "RSM DCA between positive particle and primary vertex in z direction; RSM DCA (cm); N(v0s)",100,0,1);
	fHistMCLambdacTau = new	TH1F("fHistMCLambdacTau", "Lifetime under Lambda mass hypothesis of MC lambda; Lifetime(s); N(v0s)",200,0,100);
	fHistMCLambdaNotcTau = new	TH1F("fHistMCLambdaNotcTau", "Lifetime under Lambda mass hypothesis of MC lambda background; Lifetime(s); N(v0s)",200,0,100);
	fHistMCLambdaDecayL = new	TH1F("fHistMCLambdaDecayL", "Distance between V0 and PV of MC Lambda; Distance(cm); N(v0s)",200,-0.1,30);
	fHistMCLambdaNotDecayL = new	TH1F("fHistMCLambdaNotDecayL", "Distance between V0 and PV of MC Lambda background; Distance(cm); N(v0s)",200,-0.1,30);
	fHistMcNLambdaPrimary = new	TH1F("fHistMcNLambdaPrimary","Number of primary lambdas in MC; NLambdas; i",6,-0.25,2.25);
	fHistMcNLambda = new	TH1F("fHistMcNLambda","Number of lambdas in MC; NLambdas; i",31,-0.5,30);
	fHistMcNAntilambda = new	TH1F("fHistMcNAntilambda","Number of antilambdas in MC; NAntiLambdas; i",31,-0.5,30);
	fHistMcNKshort = new	TH1F("fHistMcNKshort","Number of K0s in MC; NKshort; i",31,-0.5,30);
	fHistMK0 = new	TH1F("fHistMK0","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMLa = new	TH1F("fHistMLa","Lambda Mass; M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb = new	TH1F("fHistMLb","AntiLambda Mass; M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistNLambda = new	TH1F("fHistNLambda", "Number of lambda per event; N(lambda); N(events)",50,-0.5,49.5);
	fHistNV0 = new	TH1F("fHistNV0","V0 frequency distribution; Number of V0 Candidates",1000,0,100000);
	fHistPtV0 = new	TH1F("fHistPtV0","V0 P_{T}; P_{T} (GeV/c);dN/dP_{T} (GeV/c)^{-1}",40,0.,4.);
	fHistPVZ = new	TH1F("fHistPVZ","Z primary; Z (cm); Counts",100,-10,10);
	fHistTauLa = new	TH1F("fHistTauLa", "Lifetime under Lambda mass hypothesis; Lifetime(s); N(v0s)",200,0,100);
	fHistV0Z = new	TH1F("fHistV0Z","Z decay; Z (cm); Counts",100,-10,10);
	fHistZ = new	TH1F("fHistZ","Z decay - Z primary; Z (cm); Counts",100,-10,10);
	
	//TH2F follow
	fHistBetheBlochTPCNeg = new	TH2F("fHistBetheBlochTPCNeg","-dE/dX against Momentum for negative daughter from TPC; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);
	fHistBetheBlochTPCPos = new	TH2F("fHistBetheBlochTPCPos","-dE/dX against Momentum for positive daughter from TPC; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);
	fHistImpactxyImpactz = new	TH2F("fHistImpactxyImpactz", "RSM DCA between negative particle and primary vertex in xy plane; RSM DCA xy (cm); RSM DCA z (cm)",100,0,1,100,0,10);
	fHistMcPMK0Pt = new	TH2F("fHistMcPMK0Pt","Monte Carlo primary K0 Mass versus Pt; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPt = new	TH2F("fHistMcPMLaPt","Monte Carlo primary (& sigma0) Lambda Mass versus Pt; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPt = new	TH2F("fHistMcPMLbPt","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcV0MK0Pt = new	TH2F("fHistMcV0MK0Pt","Monte Carlo V0s passing cuts. K0 Mass versus Pt; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcV0MLaPt = new	TH2F("fHistMcV0MLaPt","Monte Carlo V0s passing cuts. Lambda Mass versus Pt; P_{perp} (GeV/c); Lambda Mass (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcV0MLbPt = new	TH2F("fHistMcV0MLbPt","Monte Carlo V0s passing cuts. Antilambda Mass versus Pt; P_{perp} (GeV/c); Antilambda Mass (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMK0Pt = new	TH2F("fHistMK0Pt","K0 Mass versus Pt; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPt = new	TH2F("fHistMLaPt","Lambda Mass versus Pt; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPt = new	TH2F("fHistMLbPt","AntiLambda Mass versus Pt; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistPtArm = new	TH2F("fHistPtArm","Podolanski-Armenteros Plot; #alpha; P_{perp} (GeV/c)",40,-1,1,80,0,0.5);
	fHistPtV0Z = new	TH2F("fHistPtV0Z","Pt of V0 vs Z position; Pt (GeV/c); Z (cm)",200,-0.1,1.9,200,-50,50);
	fHistRZ = new	TH2F("fHistRZ","R decay versus Z decay; Z (cm); R (cm)",100,-50,50,120,0,220);
	fHistXZ = new	TH2F("fHistXZ","X decay versus Z decay; Z (cm); X (cm)",100,-50,50,200,-200,200);
	fHistYZ = new	TH2F("fHistYZ","Y decay versus Z decay; Z (cm); Y (cm)",100,-50,50,200,-200,200);  
	
	
        
	// All histograms must be added to output list
        
	fOutputList->Add(fHistCosPA); 
	fOutputList->Add(fHistDCAV0Daughters); 
	fOutputList->Add(fHistDecayL); 
	fOutputList->Add(fHistImpactxyN); 
	fOutputList->Add(fHistImpactzN); 
	fOutputList->Add(fHistImpactxyP); 
	fOutputList->Add(fHistImpactzP); 
	fOutputList->Add(fHistMCLambdacTau);
	fOutputList->Add(fHistMCLambdaNotcTau);
	fOutputList->Add(fHistMCLambdaDecayL);
	fOutputList->Add(fHistMCLambdaNotDecayL);
	fOutputList->Add(fHistMcNLambdaPrimary);
	fOutputList->Add(fHistMcNLambda);
	fOutputList->Add(fHistMcNAntilambda);
	fOutputList->Add(fHistMcNKshort);
	fOutputList->Add(fHistMK0); 
	fOutputList->Add(fHistMLa); 
	fOutputList->Add(fHistMLb); 
	fOutputList->Add(fHistNLambda); 
	fOutputList->Add(fHistNV0); 
	fOutputList->Add(fHistPtV0); 
	fOutputList->Add(fHistPVZ); 
	fOutputList->Add(fHistTauLa); 
	fOutputList->Add(fHistV0Z); 
	fOutputList->Add(fHistZ);
	fOutputList->Add(fHistBetheBlochTPCNeg); 
	fOutputList->Add(fHistBetheBlochTPCPos); 
	fOutputList->Add(fHistImpactxyImpactz); 
	fOutputList->Add(fHistMcPMK0Pt);
	fOutputList->Add(fHistMcPMLaPt);
	fOutputList->Add(fHistMcPMLbPt);
	fOutputList->Add(fHistMcV0MK0Pt);
	fOutputList->Add(fHistMcV0MLaPt);
	fOutputList->Add(fHistMcV0MLbPt);
	fOutputList->Add(fHistMK0Pt); 
	fOutputList->Add(fHistMLaPt); 
	fOutputList->Add(fHistMLbPt); 
	fOutputList->Add(fHistPtArm); 
	fOutputList->Add(fHistRZ); 
	fOutputList->Add(fHistXZ); 
	fOutputList->Add(fHistYZ);
	
	
    PostData(1, fOutputList); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskLukeV0::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
	
	// paramaters used for most cuts, to minimise editing
	double cutCosPa(0.998), cutcTau(2);
	double cutNImpact(-999), cutDCA(0.4);
	double cutBetheBloch(3);
	double cutMinNClustersTPC(70), cutMaxChi2PerClusterTPC(-999);
	double isMonteCarlo(true);
	double cutEta(0.8);
	
	//Track Cuts set here
	if(cutMinNClustersTPC != -999)
	{(fTrackCuts->SetMinNClustersTPC(int(cutMinNClustersTPC)));}
	if(cutMaxChi2PerClusterTPC != -999)
	{fTrackCuts->SetMaxChi2PerClusterTPC(cutMaxChi2PerClusterTPC);}
	fTrackCuts->SetAcceptKinkDaughters(kFALSE); 
	fTrackCuts->SetRequireTPCRefit(kTRUE);
	
        
    // Create pointer to reconstructed event

	AliVEvent *event = InputEvent();
	if (!event) { Printf("ERROR: Could not retrieve event"); return; }
	
	// create pointer to event
	AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(event);
	if (!fESD) {
		AliError("Cannot get the ESD event");
		return;
	}

	/*********************************************************************/
	// MONTE CARLO SECTION
	// This section loops over all MC tracks

	int nLambdaMC = 0;
	int nAntilambdaMC = 0;
	int nKshortMC = 0;
	
	if(isMonteCarlo) 
	{

		// If the task accesses MC info, this can be done as in the commented block below:
		
		// Create pointer to reconstructed event
		AliMCEvent *mcEvent = MCEvent();
		if (!mcEvent) 
			{ 
				Printf("ERROR: Could not retrieve MC event"); 
				//return; 
			}
		else
			{
				Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
			}
			 
		// set up a stack for use in check for primary/stable particles
		AliStack* mcStack = mcEvent->Stack();
		if( !mcStack ) { Printf( "Stack not available"); return; }
		
		AliMCVertex *mcpv = (AliMCVertex *) mcEvent->GetPrimaryVertex();
		Double_t mcpvPos[3];
		if (mcpv != 0)
		{
			mcpv->GetXYZ(mcpvPos);
		}
		else 
		{
			Printf("ERROR: Could not resolve MC primary vertex");
			return;
		}
		
		//loop over all MC tracks
		for(Int_t iMCtrack = 0; iMCtrack < mcEvent->GetNumberOfTracks(); iMCtrack++)
		{
			
			//booleans to check if track is La, Lb, K0 and primary
			bool lambdaMC = false;
			bool antilambdaMC = false;
			bool kshortMC = false;
			bool isprimaryMC = false;
			
			AliMCParticle *mcPart = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iMCtrack));
			
			if(mcPart->PdgCode() == kLambda0)
			{
				lambdaMC = true;
				nLambdaMC++;
			}
			else if(mcPart->PdgCode() == kK0Short)
			{
				kshortMC = true;
				nKshortMC++;
			}
			else if(mcPart->PdgCode() == kLambda0Bar)
			{
				antilambdaMC = true;
				nAntilambdaMC++;
			}
			
			
			Int_t motherLabel = mcPart->GetMother();
			AliMCParticle *mcMother = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(motherLabel));
			double motherType = -1;
			if(motherLabel >= 0)
			{motherType = mcMother->PdgCode();}
			
			// this block of code is used to include primary Sigma0 decays as primary lambda/antilambda
			bool sigma0MC = false;
			if(motherType == 3212 || motherType == -3212)
				{
				if(mcEvent->IsPhysicalPrimary(motherLabel))
				   {sigma0MC = true;}
				}
				   
			
			if(mcEvent->IsPhysicalPrimary(iMCtrack) || sigma0MC)
			   {
				   isprimaryMC = true;
				   if(lambdaMC)
				   {
					   fHistMcNLambdaPrimary->Fill(1);
					   if(TMath::Abs(mcPart->Y())<=0.5)
						{fHistMcPMLaPt->Fill(mcPart->Pt(),mcPart->M());}
				   }
				   if(antilambdaMC)
				   {
					   if(TMath::Abs(mcPart->Y())<=0.5)
					   {fHistMcPMLbPt->Fill(mcPart->Pt(),mcPart->M());}
				   }
				   if(kshortMC)
				   {
					   if(TMath::Abs(mcPart->Y())<=0.5)
					   {fHistMcPMK0Pt->Fill(mcPart->Pt(),mcPart->M());}
				   }
				}
			else 
				{
					isprimaryMC = false;
					if(lambdaMC)
					{
						fHistMcNLambdaPrimary->Fill(2);
					}
				}
			
		}
		

	} 


	fHistMcNLambda->Fill(nLambdaMC);
	fHistMcNAntilambda->Fill(nAntilambdaMC);
	fHistMcNKshort->Fill(nKshortMC);
	
	//END OF MONTE CARLO SECTION
	/*********************************************************************/
	
			
	

	
    // Do some fast cuts first
    // check for good reconstructed vertex
    if(!(fESD->GetPrimaryVertex()->GetStatus())) return;
    // if vertex is from spd vertexZ, require more stringent cut
    if (fESD->GetPrimaryVertex()->IsFromVertexerZ()) {
        if (fESD->GetPrimaryVertex()->GetDispersion()>0.02 ||  fESD->GetPrimaryVertex()->GetZRes()>0.25 ) return; // bad vertex from VertexerZ
    }
     
	Double_t tV0Position[3];
	Double_t tPVPosition[3];
	Double_t radius;
	
	// physics selection
	UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	if(!maskIsSelected)
    {
		//printf("Event failed physics selection\n");
		return;
    }
	
	//if additional initial cuts wanted, can set conditions here
	bool isCut = (fESD->GetNumberOfTracks()==0);
	if (isCut)
	{return;}
	
	//gets primary vertex for the event
	const AliESDVertex *kPv = ((AliESDEvent *)fESD)->GetPrimaryVertex();
	if ( kPv != 0 ) 
    {
		tPVPosition[0] = kPv->GetXv();
		tPVPosition[1] = kPv->GetYv();
		tPVPosition[2] = kPv->GetZv();
		if( tPVPosition[2] == 0. ) 
		{
			//printf("WARNING: Primary vertex a Z = 0, aborting\n");
			return;
		}
    }
	else 
    {
		//printf("ERROR: Primary vertex not found\n");
		return;
    }
	if( !kPv->GetStatus())
    {return;}
	
	

	
	int nLambda(0);
	int nV0s(0);
	
	// V0 loop - runs over every v0 in the event
	for (Int_t iV0 = 0; iV0 < fESD->GetNumberOfV0s(); iV0++) 
    {
		
		AliESDv0 *v0 = fESD->GetV0(iV0);
		if (!v0) 
		{
			printf("ERROR: Could not receive v0 %d\n", iV0);
			continue;
		}
		
		bool lambdaCandidate = true;
		bool antilambdaCandidate = true;
		bool kshortCandidate = true;
		
		// keep only events of interest for fHistMLa plots
		
        if (v0->GetEffMass(4,2) < 1.08 || v0->GetEffMass(4,2) > 1.2 || TMath::Abs(v0->Y(3122))>0.5 )
		{lambdaCandidate = false;}
        if (v0->GetEffMass(2,4) < 1.08 || v0->GetEffMass(2,4) > 1.2 || TMath::Abs(v0->Y(-3122))>0.5)
		{antilambdaCandidate = false;}
        if (v0->GetEffMass(2,2) < 0.414 || v0->GetEffMass(2,2) > 0.582 || TMath::Abs(v0->Y(310))>0.5)
		{kshortCandidate = false;}
		if (v0->GetOnFlyStatus())
		{continue;}
		
		if(!isMonteCarlo) 
		{if(lambdaCandidate == false && antilambdaCandidate == false && kshortCandidate == false)
		{continue;}}
		
		
		//gets details of the v0
		v0->GetXYZ(tV0Position[0],tV0Position[1],tV0Position[2]);
		radius = TMath::Sqrt(tV0Position[0]*tV0Position[0]+tV0Position[1]*tV0Position[1]);
		
		double decayLength = (sqrt((tV0Position[0]-tPVPosition[0])*(tV0Position[0]-tPVPosition[0])+(tV0Position[1]-tPVPosition[1])*(tV0Position[1]-tPVPosition[1])+(tV0Position[2]-tPVPosition[2])*(tV0Position[2]-tPVPosition[2])));
		double cTauLa = decayLength*(v0->GetEffMass(4,2))/(v0->P());
		double cTauLb = decayLength*(v0->GetEffMass(2,4))/(v0->P());
		double cTauK0 = decayLength*(v0->GetEffMass(2,2))/(v0->P());
		
		Int_t indexP, indexN;
		indexP = TMath::Abs(v0->GetPindex());
		AliESDtrack *posTrack = ((AliESDEvent*)fESD)->GetTrack(indexP);
		indexN = TMath::Abs(v0->GetNindex());
		AliESDtrack *negTrack = ((AliESDEvent*)fESD)->GetTrack(indexN);
		
		if(!posTrack || !negTrack)
		{continue;}
		
		double pTrackMomentum[3];
		double nTrackMomentum[3];
		double pV0x, pV0y, pV0z;
		posTrack->GetConstrainedPxPyPz(pTrackMomentum);
		negTrack->GetConstrainedPxPyPz(nTrackMomentum);
		v0->GetPxPyPz(pV0x, pV0y, pV0z);

		//const double kMLambda = 1.115;
		const double kMProton = 0.938;
		//const double kMPi     = 0.140;
		
		double pPos2 = sqrt(pTrackMomentum[0]*pTrackMomentum[0]+pTrackMomentum[1]*pTrackMomentum[1]+pTrackMomentum[2]*pTrackMomentum[2]);
		double pNeg2 = sqrt(nTrackMomentum[0]*nTrackMomentum[0]+nTrackMomentum[1]*nTrackMomentum[1]+nTrackMomentum[2]*nTrackMomentum[2]);
		//double pV02 = sqrt(pV0x*pV0x+pV0y*pV0y+pV0z*pV0z);
		
		//to prevent segfaults when ratios etc taken
		//if(pV02 < 0.01 || pPos2 <0.01 || pNeg2 <0.01)
		//{continue;}
		
		Float_t pImpactxy(0), pImpactz(0);
		Float_t nImpactxy(0), nImpactz(0);
		posTrack->GetImpactParameters(pImpactxy,pImpactz);
		negTrack->GetImpactParameters(nImpactxy,nImpactz);
		nImpactxy = sqrt((nImpactxy*nImpactxy));
		nImpactz  = sqrt((nImpactz *nImpactz ));
		pImpactxy = sqrt((pImpactxy*pImpactxy));
		pImpactz  = sqrt((pImpactz *pImpactz ));
		
		/*********************************************************************/
		// Cuts are implemented here.
		
		if(!(fTrackCuts->IsSelected(posTrack)) || !(fTrackCuts->IsSelected(negTrack)))
		{
			lambdaCandidate = false;
			antilambdaCandidate = false;
			kshortCandidate = false;
		}
		
		//extra cut to account for difference between p2 & p1 data
		if(nImpactxy < 0.1 || pImpactxy < 0.1)
		{
			lambdaCandidate = false;
			antilambdaCandidate = false;
			kshortCandidate = false;
		}
		
		//psuedorapidity cut
		if(cutEta != -999)
		{
			if(TMath::Abs(posTrack->Eta()) > cutEta || TMath::Abs(negTrack->Eta())  >cutEta)
			{
				lambdaCandidate = false;
				antilambdaCandidate = false;
				kshortCandidate = false;
			}
		}
		
		//pointing angle cut
		if(cutCosPa != -999)
		{
			if (v0->GetV0CosineOfPointingAngle(tPVPosition[0],tPVPosition[1],tPVPosition[2]) < cutCosPa)
			{
				lambdaCandidate = false;
				antilambdaCandidate = false;
				kshortCandidate = false;
			}
		}
		
		//lifetime cut
		if(cutcTau != -999)
		{
			if(cTauLa < cutcTau)
			{
				lambdaCandidate = false;
			}
			if(cTauLb < cutcTau)
			{
				antilambdaCandidate = false;
			}
			if(cTauK0 < cutcTau)
			{
				kshortCandidate = false;
			}
			
		}
		
		// Impact paramater cut (on neg particle)
		if(cutNImpact != -999)
		{
			if(nImpactxy < cutNImpact || nImpactz < cutNImpact)
			{
				lambdaCandidate = false;
			}
			if(pImpactxy < cutNImpact || pImpactz < cutNImpact)
			{
				antilambdaCandidate = false;
			}
		}
		

		// DCA between daughterscut
		if(cutDCA != -999)
		{
			if(v0->GetDcaV0Daughters() > cutDCA)
			{
				lambdaCandidate = false;
				antilambdaCandidate = false;
				kshortCandidate = false;
			}
		}
		
		// Bethe Bloch cut. Made sightly complicated as options for crude cuts still included. Should probably reduce to just 'official' cuts
		if(cutBetheBloch != -999)
		{ 
			if(posTrack->GetTPCsignal() <0 || negTrack->GetTPCsignal()<0)
			{continue;}
			
			if(lambdaCandidate)
			{
				if(cutBetheBloch > 0)
				{
				if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kProton)) > cutBetheBloch )
				{lambdaCandidate = false;}
				if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kPion)) > cutBetheBloch )
				{lambdaCandidate = false;}
				}
				
				if(cutBetheBloch == -4)  
				{				
				if(isMonteCarlo) 
					{
					double beta2 = TMath::Power((pPos2/TMath::Sqrt((pPos2*pPos2+0.9*0.9))),2);
					double gamma2 = 1.0/(1.0-beta2);
					if(posTrack->GetTPCsignal() < (2.0/beta2)*(TMath::Log(1e6*beta2*gamma2)-beta2))
					{lambdaCandidate = false;}
					}
				else 
					{ 
					double beta2 = TMath::Power((pPos2/TMath::Sqrt((pPos2*pPos2+kMProton*kMProton))),2);
					double gamma2 = 1.0/(1.0-beta2);
					if(posTrack->GetTPCsignal() < (2.3/beta2)*(TMath::Log(1e6*beta2*gamma2)-beta2))
					{lambdaCandidate = false;}
					}
				}
				
			}
			
			if(antilambdaCandidate)
			{
				if(cutBetheBloch > 0)
				{
					if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kProton)) > cutBetheBloch )
					{antilambdaCandidate = false;}
					if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kPion)) > cutBetheBloch )
					{antilambdaCandidate = false;}
				}
				
				if(cutBetheBloch == -4)  
				{ 
					if(isMonteCarlo) 
					{
						double beta2 = TMath::Power((pNeg2/TMath::Sqrt((pNeg2*pNeg2+0.9*0.9))),2);
						double gamma2 = 1.0/(1.0-beta2);
						if(negTrack->GetTPCsignal() < (2.0/beta2)*(TMath::Log(1e6*beta2*gamma2)-beta2))
						{antilambdaCandidate = false;}
					}
					else 
					{ 
						double beta2 = TMath::Power((pNeg2/TMath::Sqrt((pNeg2*pNeg2+0.9*0.9))),2);
						double gamma2 = 1.0/(1.0-beta2);
						if(negTrack->GetTPCsignal() < (2.3/beta2)*(TMath::Log(1e6*beta2*gamma2)-beta2))
						{antilambdaCandidate = false;}
					}
				}
				
			}
			
			if(kshortCandidate)
			{
				if(cutBetheBloch > 0)
				{
					if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kPion)) > cutBetheBloch )
					{kshortCandidate = false;}
					if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kPion)) > cutBetheBloch )
					{kshortCandidate = false;}
				}
				
				
				if(cutBetheBloch == -4)  
				{ 
					double par0 = 0.20;
					double par1 = 4.2;
					double par2 = 1000000;
					
					if(isMonteCarlo) 
					{ 
						double beta2 = TMath::Power((pNeg2/TMath::Sqrt((pNeg2*pNeg2+par0*par0))),2);
						double gamma2 = 1.0/(1.0-beta2);
						if(negTrack->GetTPCsignal() > (par1/beta2)*(TMath::Log(par2*beta2*gamma2)-beta2) && TMath::Log10(pNeg2) > -0.6)
						{kshortCandidate = false;}
						
						beta2 = TMath::Power((pPos2/TMath::Sqrt((pPos2*pPos2+par0*par0))),2);
						gamma2 = 1.0/(1.0-beta2);
						if(posTrack->GetTPCsignal() > (par1/beta2)*(TMath::Log(par2*beta2*gamma2)-beta2) && TMath::Log10(pNeg2) > -0.6)
						{kshortCandidate = false;}
					}
					else 
					{ 
						double beta2 = TMath::Power((pNeg2/TMath::Sqrt((pNeg2*pNeg2+par0*par0))),2);
						double gamma2 = 1.0/(1.0-beta2);
						if(negTrack->GetTPCsignal() > (par1/beta2)*(TMath::Log(par2*beta2*gamma2)-beta2) && TMath::Log10(pNeg2) > -0.6)
						{kshortCandidate = false;}
						
						beta2 = TMath::Power((pPos2/TMath::Sqrt((pPos2*pPos2+par0*par0))),2);
						gamma2 = 1.0/(1.0-beta2);
						if(posTrack->GetTPCsignal() > (par1/beta2)*(TMath::Log(par2*beta2*gamma2)-beta2) && TMath::Log10(pPos2) > -0.6)
						{kshortCandidate = false;}
					}
				}
				
			}
		}
		
		// Selection of random cuts which I've been playing with
		/*if(nImpactxy > 3 || pImpactxy > 2)
		 {				
		 lambdaCandidate = false;
		 }*/
		
		/*if(nImpactxy < 0.4 || pImpactxy < 0.4 || nImpactxy > 2.5 || pImpactxy >2.5)
		 {				
		 antilambdaCandidate = false;
		 }	*/
		
		if(decayLength > 15 )
		{lambdaCandidate = false;}
		
		// Cuts to look at just lowpt region of lambdas
		//if(v0->Pt() < 0.3 || v0->Pt() > 0.7 || !lambdaCandidate)
		//{continue;}
		
		// cuts to just look at signal/background region of lambda mass
		//if(!((v0->GetEffMass(4,2) >= 1.096 && v0->GetEffMass(4,2) < 1.106) || (v0->GetEffMass(4,2) >= 1.126 && v0->GetEffMass(4,2) < 1.136) ))
		//if(!(v0->GetEffMass(4,2) >= 1.106 && v0->GetEffMass(4,2) < 1.126 ))
		//{continue;}
		/*********************************************************************/

		/*********************************************************************/
		// MONTE CARLO SECTION 2
		// this section looks at the individual V0s
		
		bool mcLambdaCandidate = true;
		bool mcAntilambdaCandidate = true;
		bool mcK0Candidate = true;
		bool realParticle = true;
		
		if(isMonteCarlo) 
		{
			
			AliMCEvent *mcEvent = MCEvent();
			AliStack* mcStack = mcEvent->Stack();
			if( !mcStack ) { Printf( "Stack not available"); return; }
			
			TParticle *negParticle = mcStack->Particle( TMath::Abs(negTrack->GetLabel())); 
			TParticle *posParticle = mcStack->Particle( TMath::Abs(posTrack->GetLabel())); 
			
			Int_t negParticleMotherLabel = negParticle->GetFirstMother();
			Int_t posParticleMotherLabel = posParticle->GetFirstMother();
			
			if( negParticleMotherLabel == -1 || posParticleMotherLabel == -1)
			{
			realParticle = false;
			mcLambdaCandidate = false;
			mcAntilambdaCandidate = false;
			mcK0Candidate  =false;
			}

			if( negParticleMotherLabel != posParticleMotherLabel)
			{
				mcLambdaCandidate = false;
				mcAntilambdaCandidate = false;
				mcK0Candidate  =false;
			}
			
			if(realParticle == true)
			{
				AliMCParticle *mcPart2 = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(negParticleMotherLabel));

				if(mcPart2->PdgCode() != kLambda0)
					{mcLambdaCandidate = false;}
				if(mcPart2->PdgCode() != kLambda0Bar)
					{mcAntilambdaCandidate = false;}
				if(mcPart2->PdgCode() != kK0Short)
					{mcK0Candidate  =false;}

				if(mcLambdaCandidate && lambdaCandidate)
				{
					fHistMcV0MLaPt->Fill(v0->Pt(),v0->GetEffMass(4,2)); 
				}
				if(mcAntilambdaCandidate && antilambdaCandidate)
				{
					fHistMcV0MLbPt->Fill(v0->Pt(),v0->GetEffMass(2,4));
				}
				if(mcK0Candidate && kshortCandidate)
				{
					fHistMcV0MK0Pt->Fill(v0->Pt(),v0->GetEffMass(2,2));
				}
				}

			if(mcLambdaCandidate && lambdaCandidate)
			{
				fHistMCLambdacTau->Fill(cTauLa);
				fHistMCLambdaDecayL->Fill(decayLength);
			}
			if(!mcLambdaCandidate && lambdaCandidate)
			{
			fHistMCLambdaNotcTau->Fill(cTauLa);
			fHistMCLambdaNotDecayL->Fill(decayLength);
			}
					
			}
		
	
		// END OF MONTE CARLO SECTION 2
		/*********************************************************************/

		//remove all non-candidates
		if(lambdaCandidate == false && antilambdaCandidate == false && kshortCandidate == false)
		{continue;}
		
		
		//count number of valid v0s
		nV0s+=1;
			
		/*********************************************************************/
		//This section fills histograms
		
		fHistDCAV0Daughters->Fill(v0->GetDcaV0Daughters());
		fHistCosPA->Fill(v0->GetV0CosineOfPointingAngle(tPVPosition[0],tPVPosition[1],tPVPosition[2]));
		fHistDecayL->Fill(decayLength);
		fHistTauLa->Fill(cTauLa);
		
		fHistBetheBlochTPCPos->Fill(TMath::Log10(pPos2),posTrack->GetTPCsignal());
		fHistBetheBlochTPCNeg->Fill(TMath::Log10(pNeg2),negTrack->GetTPCsignal());

		fHistImpactxyN->Fill(nImpactxy);
		fHistImpactzN->Fill(nImpactz);
		fHistImpactxyP->Fill(pImpactxy);
		fHistImpactzP->Fill(pImpactz);
		
		fHistImpactxyImpactz->Fill(nImpactxy,nImpactz);
		
		fHistV0Z->Fill(tV0Position[2]);
		fHistZ->Fill(tV0Position[2]-tPVPosition[2]);
	
		fHistRZ->Fill(tV0Position[2],radius);
		fHistPtV0Z->Fill(v0->Pt(),tV0Position[2]);
		
		fHistPtArm->Fill(v0->AlphaV0(),v0->PtArmV0());
		fHistXZ->Fill(tV0Position[2],tV0Position[0]);
		fHistYZ->Fill(tV0Position[2],tV0Position[1]);
		fHistPtV0->Fill(v0->Pt());
		
		//effective mass histograms
		
		//sets assumed particle type of pos/neg daughters. 
		// 0 = electron, 1 = Muon, 2 = pion, 3 = kaon, 4 = proton.
		int dPos = 0;
		int dNeg = 0;
		
		//    v0->ChangeMassHypothesis(kLambda0);
		dPos = 4;
		dNeg = 2;
		if(v0->GetEffMass(dPos,dNeg) > 1.11 && v0->GetEffMass(dPos,dNeg) < 1.13)
		{
			if(!(v0->GetOnFlyStatus()))
			{
				nLambda++;
			}
		}
		if(lambdaCandidate)
		{
			fHistMLa->Fill(v0->GetEffMass(dPos,dNeg));  
			fHistMLaPt->Fill(v0->Pt(),v0->GetEffMass(dPos,dNeg));   
		}
		
		//    v0->ChangeMassHypothesis(kK0Short);
		dPos = 2;
		dNeg = 2;
		if(kshortCandidate)
		{
			fHistMK0->Fill(v0->GetEffMass(dPos,dNeg));
			fHistMK0Pt->Fill(v0->Pt(),v0->GetEffMass(dPos,dNeg));
		}
		//    v0->ChangeMassHypothesis(kLambda0Bar);
		dPos = 2;
		dNeg = 4;
		if(antilambdaCandidate)
		{
			fHistMLb->Fill(v0->GetEffMass(dPos,dNeg));
			fHistMLbPt->Fill(v0->Pt(),v0->GetEffMass(dPos,dNeg));
		}
				
    }
	
	
	fHistPVZ->Fill(tPVPosition[2]);
	fHistNV0->Fill(nV0s);
	fHistNLambda->Fill(nLambda);	
	
	// NEW HISTO should be filled before this point, as PostData puts the
    // information for this iteration of the UserExec in the container
	PostData(1, fOutputList);
}


//________________________________________________________________________
void AliAnalysisTaskLukeV0::Terminate(Option_t *) 
{
    // Draw result to screen, or perform fitting, normalizations
    // Called once at the end of the query
        
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if(!fOutputList) { Printf("ERROR: could not retrieve TList fOutputList"); return; }
      
 }
