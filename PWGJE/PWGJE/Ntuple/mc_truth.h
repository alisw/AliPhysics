#include <AliMCEvent.h>

namespace {

    bool final_state_primary(AliMCEvent *mc_event, Int_t index)
    {
        if (mc_event->HasSubsidiaries()) {
            AliMCEvent *e = NULL;
            Int_t j = mc_event->FindIndexAndEvent(index, e);
            return final_state_primary(e, j);
        }
        else if (mc_event != NULL) {
            AliStack *s = mc_event->Stack();

            return index < s->GetNprimary() &&
                s->Particle(index)->GetStatusCode() == 1;
        }
        else {
            return false;
        }
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

    bool parton_cms_algorithmic(AliMCEvent *mc_event, Int_t index)
    {
        AliStack *s = mc_event->Stack();
        const AliMCParticle *p =
            dynamic_cast<AliMCParticle *>(mc_event->GetTrack(index));

        if (!(s != NULL && p != NULL)) {
            return false;
        }

        return pdg_is_parton(p->PdgCode()) &&
            !(p->GetFirstDaughter() > 0 &&
              p->GetFirstDaughter() < s->GetNprimary() &&
              dynamic_cast<AliMCParticle *>(
                mc_event->GetTrack(
                    p->GetFirstDaughter())) != NULL &&
              pdg_is_parton(dynamic_cast<AliMCParticle *>(
                mc_event->GetTrack(
                    p->GetFirstDaughter()))->PdgCode())) &&
            !(p->GetLastDaughter() > 0 &&
              p->GetLastDaughter() < s->GetNprimary() &&
              dynamic_cast<AliMCParticle *>(
                mc_event->GetTrack(
                    p->GetLastDaughter())) != NULL &&
              pdg_is_parton(dynamic_cast<AliMCParticle *>(
                mc_event->GetTrack(
                    p->GetLastDaughter()))->PdgCode())) &&
            !(p->Px() == 0 && p->Py() == 0);
    }

}
