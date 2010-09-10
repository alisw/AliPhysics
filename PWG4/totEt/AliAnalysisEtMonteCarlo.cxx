#include "AliAnalysisEtMonteCarlo.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TH2F.h"

Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliVEvent* ev)
{
     ResetEventValues();
     
    // Get us an mc event
    AliMCEvent *event = dynamic_cast<AliMCEvent*>(ev);

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

        // Check for reasonable (for now neutral and singly charged) charge on the particle
        //TODO:Maybe not only singly charged?
        if (TMath::Abs(pc->Charge()) != EtMonteCarloCuts::kSingleChargedParticle && pc->Charge() != EtMonteCarloCuts::kNeutralParticle) continue;

        fMultiplicity++;

        if (TMath::Abs(part->Eta()) < fEtaCut)
        {

	  TParticlePDG *pdgCode =  part->GetPDG(0);
            if (
                TMath::Abs(pdgCode->PdgCode()) == ProtonCode ||
                TMath::Abs(pdgCode->PdgCode()) == NeutronCode ||
                TMath::Abs(pdgCode->PdgCode()) == LambdaCode ||
                TMath::Abs(pdgCode->PdgCode()) == XiCode ||
                TMath::Abs(pdgCode->PdgCode()) == Xi0Code ||
                TMath::Abs(pdgCode->PdgCode()) == OmegaCode
	       )
            {
                particleMassPart = -TMath::Sign(pdgCode->PdgCode(), pdgCode->PdgCode())*pdgCode->Mass();
            }
            
            if (pdgCode->Charge() == EtMonteCarloCuts::kNeutralParticle)
            {
	       fNeutralMultiplicity++;
                fTotNeutralEt += part->Energy()*TMath::Sin(part->Theta());

                if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
                {
                    fTotNeutralEtAcc += part->Energy()*TMath::Sin(part->Theta());
                    fTotEtAcc += part->Energy()*TMath::Sin(part->Theta());
                }
            }
            else if (pdgCode->Charge() != EtMonteCarloCuts::kNeutralParticle)
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
{

    AliAnalysisEt::Init();

    fVertexXCut = EtReconstructedCuts::kVertexXCut;
    fVertexYCut = EtReconstructedCuts::kVertexYCut;
    fVertexZCut = EtReconstructedCuts::kVertexZCut;
    fIPxyCut = EtReconstructedCuts::kIPxyCut;
    fIPzCut = EtReconstructedCuts::kIPzCut;

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

