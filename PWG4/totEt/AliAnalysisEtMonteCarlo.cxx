//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//  implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEtMonteCarlo.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TH2F.h"
#include "TParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"

using namespace std;

ClassImp(AliAnalysisEtMonteCarlo);


// ctor
AliAnalysisEtMonteCarlo::AliAnalysisEtMonteCarlo() :
  AliAnalysisEt()
  ,fImpactParameter(0)
  ,fNcoll(0)
  ,fNpart(0)
{
}

// dtor
AliAnalysisEtMonteCarlo::~AliAnalysisEtMonteCarlo() 
{
}

Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliVEvent* ev)
{ // analyse MC event
     ResetEventValues();
     
    // Get us an mc event
     if(!ev){
            Printf("ERROR: Event does not exist");   
	    return 0;
     }
     AliMCEvent *event = dynamic_cast<AliMCEvent*>(ev);

    Double_t protonMass =fgProtonMass;

    // Hijing header
    AliGenEventHeader* genHeader = event->GenEventHeader();
    if(!genHeader){
            Printf("ERROR: Event generation header does not exist");   
	    return 0;
    }
    AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
    if (hijingGenHeader) {
      fImpactParameter = hijingGenHeader->ImpactParameter();
      fNcoll = hijingGenHeader->HardScatters(); // or should this be some combination of NN() NNw() NwN() NwNw() ?
      fNpart = hijingGenHeader->ProjectileParticipants() + hijingGenHeader->TargetParticipants(); 
      /*
      printf("Hijing: ImpactParameter %g ReactionPlaneAngle %g \n",
	     hijingGenHeader->ImpactParameter(), hijingGenHeader->ReactionPlaneAngle());
      printf("HardScatters %d ProjecileParticipants %d TargetParticipants %d\n",
	     hijingGenHeader->HardScatters(), hijingGenHeader->ProjectileParticipants(), hijingGenHeader->TargetParticipants()); 
      printf("ProjSpectatorsn %d ProjSpectatorsp %d TargSpectatorsn %d TargSpectatorsp %d\n",
	     hijingGenHeader->ProjSpectatorsn(), hijingGenHeader->ProjSpectatorsp(), hijingGenHeader->TargSpectatorsn(), hijingGenHeader->TargSpectatorsp());
      printf("NN %d NNw %d NwN %d, NwNw %d\n",
	     hijingGenHeader->NN(), hijingGenHeader->NNw(), hijingGenHeader->NwN(), hijingGenHeader->NwNw());
      */
    }

    /* // placeholder if we want to get some Pythia info later
    AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
    if (pythiaGenHeader) { // not Hijing; try with Pythia      
      printf("Pythia: ProcessType %d  GetPtHard %g \n",
	     pythiaGenHeader->ProcessType(), pythiaGenHeader->GetPtHard());
    }
    */

    // Let's play with the stack!
    AliStack *stack = event->Stack();

    Int_t nPrim = stack->GetNtrack();

    for (Int_t iPart = 0; iPart < nPrim; iPart++)
    {

        TParticle *part = stack->Particle(iPart);

        if (!part)
        {
            Printf("ERROR: Could not get particle %d", iPart);
            continue;
        }

        TParticlePDG *pdg = part->GetPDG(0);

	Double_t particleMassPart = 0; //The mass part in the Et calculation for this particle

        // Check if it is a primary particle
        if (!stack->IsPhysicalPrimary(iPart)) continue;

	// printf("MC: iPart %03d eta %4.3f phi %4.3f code %d charge %g \n", iPart, part->Eta(), part->Phi(), pdg->PdgCode(), pdg->Charge()); // tmp/debug printout

        // Check for reasonable (for now neutral and singly charged) charge on the particle
        //TODO:Maybe not only singly charged?
        if (TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloSingleChargedParticle())<1e-3 && TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloNeutralParticle())<1e-3) continue;

        fMultiplicity++;

        if (TMath::Abs(part->Eta()) < fCuts->GetCommonEtaCut())
        {

            if (
                TMath::Abs(pdg->PdgCode()) == fgProtonCode ||
                TMath::Abs(pdg->PdgCode()) == fgNeutronCode ||
                TMath::Abs(pdg->PdgCode()) == fgLambdaCode ||
                TMath::Abs(pdg->PdgCode()) == fgXiCode ||
                TMath::Abs(pdg->PdgCode()) == fgXi0Code ||
                TMath::Abs(pdg->PdgCode()) == fgOmegaCode
	       )
            {
	      if (pdg->PdgCode() > 0) { particleMassPart = - protonMass;}
	      if (pdg->PdgCode() < 0) { particleMassPart = protonMass;}
	    }
	    Double_t et = part->Energy() * TMath::Sin(part->Theta()) + particleMassPart;
	    	      
	    if (pdg->PdgCode() == fgProtonCode || pdg->PdgCode() == fgAntiProtonCode)
	      {
		fProtonEt += et;
	      }
	    if (pdg->PdgCode() == fgPiPlusCode || pdg->PdgCode() == fgPiMinusCode)
	      {
		fPionEt += et;
	      }
	    if (pdg->PdgCode() == fgKPlusCode || pdg->PdgCode() == fgKMinusCode)
	      {
		fChargedKaonEt += et;
	      }
	    if (pdg->PdgCode() == fgMuPlusCode || pdg->PdgCode() == fgMuMinusCode)
	      {
		fMuonEt += et;
	      }
	    if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
	      {
		fElectronEt += et;
	      }

	    // some neutrals also
	    if(pdg->PdgCode() == fgNeutronCode)
	    {
	      fNeutronEt += et;
	    }
            if(pdg->PdgCode() == fgAntiNeutronCode)
	    {
	      fAntiNeutronEt += et;
	    }
	    if(pdg->PdgCode() == fgGammaCode)
	    {
	      fGammaEt += et;
	    }

            if (TMath::Abs(pdg->Charge() - fCuts->GetMonteCarloNeutralParticle()) <1e-3 )
            {
	       fNeutralMultiplicity++;
	       fTotNeutralEt += et;

                if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
                {
		  fTotNeutralEtAcc += et;
		  fTotEtAcc += et;
                }
            }
            else if (TMath::Abs( pdg->Charge() - fCuts->GetMonteCarloNeutralParticle())>1e-3 )
            {
	       fChargedMultiplicity++;
	       fTotChargedEt += et;
                if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
                {
		  fTotChargedEtAcc += et;
		  fTotEtAcc += et;


		    if (pdg->PdgCode() == fgProtonCode || pdg->PdgCode() == fgAntiProtonCode)
		      {
			fProtonEtAcc += et;
		      }
		    if (pdg->PdgCode() == fgPiPlusCode || pdg->PdgCode() == fgPiMinusCode)
		      {
			fPionEtAcc += et;
		      }
		    if (pdg->PdgCode() == fgKPlusCode || pdg->PdgCode() == fgKMinusCode)
		      {
			fChargedKaonEtAcc += et;
		      }
		    if (pdg->PdgCode() == fgMuPlusCode || pdg->PdgCode() == fgMuMinusCode)
		      {
			fMuonEtAcc += et;
		      }
		    if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
		      {
			fElectronEtAcc += et;
		      }
		    
                }

		//	  if (TrackHitsCalorimeter(part, event->GetMagneticField()))
		if (TrackHitsCalorimeter(part)) // magnetic field info not filled?
		  {
		    if (pdg->Charge() > 0) fHistPhivsPtPos->Fill(part->Phi(),part->Pt());
		    else fHistPhivsPtNeg->Fill(part->Phi(), part->Pt());
		  }
	    }
	}
    }
    
    fTotEt = fTotChargedEt + fTotNeutralEt;
    fTotEtAcc = fTotChargedEtAcc + fTotNeutralEtAcc;
    
    FillHistograms();

    return 0;    
}

void AliAnalysisEtMonteCarlo::Init()
{ // init
    AliAnalysisEt::Init();
}

void AliAnalysisEtMonteCarlo::ResetEventValues()
{ // reset event values
  AliAnalysisEt::ResetEventValues();

  // collision geometry defaults for p+p:
  fImpactParameter = 0;
  fNcoll = 1;
  fNpart = 2;  
}

void AliAnalysisEtMonteCarlo::CreateHistograms()
{ // histogram related additions
  AliAnalysisEt::CreateHistograms();
  if (fTree) {
    fTree->Branch("fImpactParameter",&fImpactParameter,"fImpactParameter/D");
    fTree->Branch("fNcoll",&fNcoll,"fNcoll/I");
    fTree->Branch("fNpart",&fNpart,"fNpart/I");
  }
}

bool AliAnalysisEtMonteCarlo::TrackHitsCalorimeter(TParticle* part, Double_t magField)
{
  //  printf(" TrackHitsCalorimeter - magField %f\n", magField);
   AliESDtrack *esdTrack = new AliESDtrack(part);
   // Printf("MC Propagating track: eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());

    Bool_t prop = esdTrack->PropagateTo(fDetectorRadius, magField);

    // if(prop) Printf("Track propagated, eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());

    bool status = prop && 
      TMath::Abs(esdTrack->Eta()) < fEtaCutAcc && 
      esdTrack->Phi() > fPhiCutAccMin*TMath::Pi()/180. && 
      esdTrack->Phi() < fPhiCutAccMax*TMath::Pi()/180.;
    delete esdTrack;

    return status;
}

