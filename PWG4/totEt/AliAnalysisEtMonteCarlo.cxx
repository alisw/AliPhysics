#include "AliAnalysisEtMonteCarlo.h"
#include "AliAnalysisEtCuts.h"

#include "AliStack.h"
#include "AliMCEvent.h"

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
                TMath::Abs(pdgCode->PdgCode()) == fPdgDB->GetParticle("proton")->PdgCode() ||
                TMath::Abs(pdgCode->PdgCode()) == fPdgDB->GetParticle("neutron")->PdgCode() ||
                TMath::Abs(pdgCode->PdgCode()) == fPdgDB->GetParticle("Lambda0")->PdgCode() ||
                TMath::Abs(pdgCode->PdgCode()) == fPdgDB->GetParticle("Xi-")->PdgCode() ||
                TMath::Abs(pdgCode->PdgCode()) == fPdgDB->GetParticle("Xi0")->PdgCode() ||
                TMath::Abs(pdgCode->PdgCode()) == fPdgDB->GetParticle("Omega-")->PdgCode()
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
