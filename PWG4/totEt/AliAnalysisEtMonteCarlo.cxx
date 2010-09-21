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
    AliMCEvent *event = dynamic_cast<AliMCEvent*>(ev);

    // Hijing header
    AliGenEventHeader* genHeader = event->GenEventHeader();
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

    Double_t particleMassPart = 0; //The mass part in the Et calculation for this particle

    for (Int_t iPart = 0; iPart < nPrim; iPart++)
    {

        TParticle *part = stack->Particle(iPart);

        if (!part)
        {
            Printf("ERROR: Could not get particle %d", iPart);
            continue;
        }

        TParticlePDG *pc = part->GetPDG(0);

        // Check if it is a primary particle
        if (!stack->IsPhysicalPrimary(iPart)) continue;

	// printf("MC: iPart %03d eta %4.3f phi %4.3f code %d charge %g \n", iPart, part->Eta(), part->Phi(), pc->PdgCode(), pc->Charge()); // tmp/debug printout

        // Check for reasonable (for now neutral and singly charged) charge on the particle
        //TODO:Maybe not only singly charged?
        if (TMath::Abs(pc->Charge()) != fCuts->GetMonteCarloSingleChargedParticle() && pc->Charge() != fCuts->GetMonteCarloNeutralParticle()) continue;

        fMultiplicity++;

        if (TMath::Abs(part->Eta()) < fCuts->GetCommonEtaCut())
        {

	  TParticlePDG *pdgCode =  part->GetPDG(0);
            if (
                TMath::Abs(pdgCode->PdgCode()) == fProtonCode ||
                TMath::Abs(pdgCode->PdgCode()) == fNeutronCode ||
                TMath::Abs(pdgCode->PdgCode()) == fLambdaCode ||
                TMath::Abs(pdgCode->PdgCode()) == fXiCode ||
                TMath::Abs(pdgCode->PdgCode()) == fXi0Code ||
                TMath::Abs(pdgCode->PdgCode()) == fOmegaCode
	       )
            {
                particleMassPart = -TMath::Sign(pdgCode->PdgCode(), pdgCode->PdgCode())*pdgCode->Mass();
            }
            
            if (pdgCode->Charge() == fCuts->GetMonteCarloNeutralParticle() )
            {
	       fNeutralMultiplicity++;
                fTotNeutralEt += part->Energy()*TMath::Sin(part->Theta());

                if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
                {
                    fTotNeutralEtAcc += part->Energy()*TMath::Sin(part->Theta());
                    fTotEtAcc += part->Energy()*TMath::Sin(part->Theta());
                }
            }
            else if (pdgCode->Charge() != fCuts->GetMonteCarloNeutralParticle() )
            {
	       fChargedMultiplicity++;
                fTotChargedEt += part->Energy()*TMath::Sin(part->Theta());
                if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
                {
                    fTotChargedEtAcc += part->Energy()*TMath::Sin(part->Theta());
                    fTotEtAcc += part->Energy()*TMath::Sin(part->Theta());
                }

		//	  if (TrackHitsCalorimeter(part, event->GetMagneticField()))
		if (TrackHitsCalorimeter(part)) // magnetic field info not filled?
		  {
		    if (pdgCode->Charge() > 0) fHistPhivsPtPos->Fill(part->Phi(),part->Pt());
		    else fHistPhivsPtNeg->Fill(part->Phi(), part->Pt());
		  }
	    }
	}
    }
    
    fTotNeutralEtAcc = fTotNeutralEt;
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

