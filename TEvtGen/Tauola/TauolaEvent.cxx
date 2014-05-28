#include "TauolaEvent.h"
#include "Plots.h"

using namespace std;

namespace Tauolapp
{

void TauolaEvent::undecayTaus(){
  
  std::vector<TauolaParticle*> particle_list;
  particle_list = findParticles(Tauola::getDecayingParticle());

  for(int p=0; p < (int) particle_list.size(); p++)
    particle_list.at(p)->findLastSelf()->undecay();

}

void TauolaEvent::decayTaus(){

  std::vector<TauolaParticle*> particle_list;
  particle_list = findStableParticles(Tauola::getDecayingParticle());

  while(particle_list.size()!=0){

    // tau and its matching tau-like partner is removed from the list here:
    TauolaParticlePair t_pair(particle_list);

    //t_pair.print();
    t_pair.decayTauPair();
    t_pair.checkMomentumConservation();
  }

  // Final event record modifications
  eventEndgame();
}

} // namespace Tauolapp
