/*
  Ludmila Malinina  malinina@lav01.sinp.msu.ru,   SINP MSU/Moscow and JINR/Dubna
  Ionut Arsene  i.c.arsene@fys.uio.no,            Oslo University and ISS-Bucharest
  Date        : 2007/05/30
*/

#include "TVector3.h"
#ifndef INITIAL_STATE
#include "InitialState.h"
#endif
#ifndef HADRONDECAYER_INCLUDED
#include "HadronDecayer.h"
#endif
#include <iostream> 
#include <fstream>
using namespace std;
// aic(2008/08/08): define the fLastIndex static variable 
Int_t Particle::fLastIndex;

void InitialState::Evolve(List_t &source, List_t &secondaries, ParticleAllocator &allocator, Double_t weakDecayLimit) {
  cout << "InitialState::Evolve()  IN" << endl;
  LPIT_t it;
  LPIT_t e;

  // Initialize the fLastIndex to -1
  it = source.begin();
  it->InitIndexing();

  //  cout << "InitialState::Evolve()  Loop over primaries" << endl;
  // Loop over source list (primary particles) 
  for(it = source.begin(), e = source.end(); it != e; ++it) {
    //    cout << "InitialState::Evolve()  particle(in primary loop) code = " << it->Encoding() << endl;
    //    if(GetPDGParticleStableStatus(it->Encoding())) {
      //      cout << "InitialState::Evolve() its a stable marked particle, decay time should be 0" << endl;
      //    }
    // calculate the decay time of the particle
    Double_t decayTime = GetDecayTime(*it, weakDecayLimit);
    //    cout << "InitialState::Evolve()  decay time = " << decayTime << endl;
    TVector3 shift(it->Mom().BoostVector());
    shift *= decayTime;
    it->SetLastInterTime(it->T() + decayTime);

    // if the particle is unstable run the decayer 
    if(decayTime > 0.) {
      //      cout << "InitialState::Evolve()  perform decay..." << endl;
      it->Pos(it->Pos() += TLorentzVector(shift, 0.));
      it->T(it->T() + decayTime);
      // perform the decay procedure
      // the primaries are added to the list of secondaries inside Decay()
      Decay(secondaries, *it, allocator, fDatabase);  
    }
    // else just add the primary to the list of particles
    else {
      //      cout << "InitialState::Evolve()  stable particle (just add it to the list of final particles)" << endl;
      it->SetIndex();
      allocator.AddParticle(*it, secondaries);
    }
  }

  //  cout << "InitialState::Evolve()  Loop over secondaries" << endl;
  // Loop over the secondaries list and decay the cascades
  for(it = secondaries.begin(), e = secondaries.end(); it != e;) {
    //    cout << "InitialState::Evolve()  particle(in secondaries loop) code " << it->Encoding() << endl;
    if(GetPDGParticleStableStatus(it->Encoding())) {
      //      cout << "InitialState::Evolve()  its a stable particle " << it->Encoding()
      //	   << "(already added to final list)"  << endl;
      ++it;
      continue;
    }
    // the primaries are already decayed, so just ignore it
    if(it->GetMother()==-1) {  
      //      cout << "InitialState::Evolve()  its a primary (already decayed)" << endl;
      ++it;
      continue;
    }
    Double_t decayTime = GetDecayTime(*it, weakDecayLimit);
    //    cout << "InitialState::Evolve()  decay time = " << decayTime << endl;
    it->SetLastInterTime(it->T() + decayTime);
    TVector3 shift(it->Mom().BoostVector());
    shift *= decayTime; 
    
    // if the particle is unstable run the decays
    if(decayTime > 0.) {
      //      cout << "InitialState::Evolve()  perform decay..." << endl;
      it->Pos(it->Pos() += TLorentzVector(shift, 0.));
      it->T(it->T() + decayTime);
      // The decay products are added to the list inside the Decay() function
      // The mother particle is not removed from list but can be distinguished 
      // since the daughter indexes will be initialized
      Decay(secondaries, *it, allocator, fDatabase);
      ++it;
    }
    // else just continue
    else {
      //      cout << "InitialState::Evolve() no decay, continue to next particle" << endl;
      ++it;
    }
  }
  //  cout << "InitialState::Evolve()  OUT"<< endl;
}
