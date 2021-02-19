namespace {

    bool final_state_primary(AliMCParticleContainer *mc_container, Int_t index)
    {
        const AliAODMCParticle *p = mc_container->GetMCParticle(index);

        if (p == NULL) {
            return false;
        }

        return p->MCStatusCode() == 1;
    }

    bool pdg_is_parton(Int_t pdg_code)
    {
        const int pdg_code_particle_unexcited =
            std::abs(static_cast<int>(pdg_code)) % 100000;

        switch (pdg_code_particle_unexcited) {
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 21:
            return true;
            break;
        }

        return false;
    }

    // Port of the CMS b-tagging group's "algorithmic definition" of
    // partons (at the end of showering) for PYTHIA 8. See
    // https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/
    // JetMCAlgos/src/Pythia8PartonSelector.cc

    bool parton_cms_algorithmic(AliMCParticleContainer *mc_container, Int_t index)
    {
        const AliAODMCParticle *p = mc_container->GetMCParticle(index);

        if (p == NULL) {
            return false;
        }

        return pdg_is_parton(p->PdgCode()) &&
            !(p->GetDaughterFirst() > 0 &&
              p->GetDaughterFirst() < mc_container->GetNParticles() &&
                mc_container->GetMCParticle(
                    p->GetDaughterFirst()) != NULL &&
              pdg_is_parton(
                mc_container->GetMCParticle(
                    p->GetDaughterFirst())->PdgCode())) &&
            !(p->GetDaughterLast() > 0 &&
              p->GetDaughterLast() < mc_container->GetNParticles() &&
                mc_container->GetMCParticle(
                    p->GetDaughterLast()) != NULL &&
              pdg_is_parton(
                mc_container->GetMCParticle(
                    p->GetDaughterLast())->PdgCode())) &&
            !(p->Px() == 0 && p->Py() == 0);
    }

}
