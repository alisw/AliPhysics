#include "TauolaHepMCEvent.h"
#include "Log.h"

using namespace std;

namespace Tauolapp
{

TauolaHepMCEvent::TauolaHepMCEvent(HepMC::GenEvent * event){
  m_event=event;

  // Default units
  m_momentum_unit = "GEV";
  m_length_unit   = "MM";

  if(m_event->momentum_unit() != HepMC::Units::GEV) m_momentum_unit = "MEV";
  if(m_event->length_unit()   != HepMC::Units::MM ) m_length_unit   = "CM";

  // If needed - change units used by HepMC to GEV and MM
  if( m_event->momentum_unit() != HepMC::Units::GEV ||
      m_event->length_unit()   != HepMC::Units::MM     )
  {
    m_event->use_units(HepMC::Units::GEV,HepMC::Units::MM);
  }
}

TauolaHepMCEvent::~TauolaHepMCEvent(){

  while(m_tau_list.size()!=0){
    TauolaParticle * temp = m_tau_list.back();
    m_tau_list.pop_back();
    delete temp;
  }

}

HepMC::GenEvent * TauolaHepMCEvent::getEvent(){
  return m_event;
}

std::vector<TauolaParticle*> TauolaHepMCEvent::findParticles(int pdg_id){
  
  if(m_tau_list.size()==0){

    HepMC::GenEvent::particle_const_iterator part_itr = m_event->particles_begin();
    //loop over all particle in the event looking for taus (or other)
    for( ; part_itr!=m_event->particles_end(); part_itr++){
      if(abs((*part_itr)->pdg_id())==pdg_id)
        m_tau_list.push_back(new TauolaHepMCParticle(*part_itr));
    }
  }
  return m_tau_list;
}

std::vector<TauolaParticle*> TauolaHepMCEvent::findStableParticles(int pdg_id){

  /**  HepMC::GenEvent::particle_const_iterator part_itr = m_event->particles_begin();
  //loop over all particle in the event looking for taus (or other)
  for( ; part_itr!=m_event->particles_end(); part_itr++){
    if(fabs((*part_itr)->pdg_id())==pdg_id){
      if((*part_itr)->end_vertex()){
        cout << "WARNING: Particle with pdg code " << (*part_itr)->pdg_id()
             << " has end vertex" <<endl;
      }
      else
        list.push_back(new TauolaHepMCParticle(*part_itr));
    }
    }**/

  std::vector<TauolaParticle*> tau_list = findParticles(pdg_id);
  std::vector<TauolaParticle*> stable_tau_list;

  for(int i=0; i<(int) tau_list.size(); i++){
    
    if(!tau_list.at(i)->hasDaughters())
      stable_tau_list.push_back(tau_list.at(i)); 
    else
    {
      std::vector<TauolaParticle*> t = tau_list.at(i)->getDaughters();
      //Ignore taus that we won't be decaying anyway
      if(t.size()==1) continue;
      if(t.size()==2 && (abs(t[0]->getPdgID())==15 || abs(t[1]->getPdgID())==15) ) continue;
      Log::Warning()<<"Particle with pdg code "<<tau_list.at(i)->getPdgID()
                    <<" already has daughters" <<endl;
    }
  }

  return stable_tau_list;

}

void TauolaHepMCEvent::eventEndgame(){

  //Set output units for the event
  string momentum("GEV"),length("MM");

  switch(Tauola::momentumUnit)
  {
    case Tauola::GEV:
      momentum = "GEV";
      break;
    case Tauola::MEV:
      momentum = "MEV";
      break;
    default:
      momentum = m_momentum_unit;
  }

  switch(Tauola::lengthUnit)
  {
    case Tauola::MM:
      length = "MM";
      break;
    case Tauola::CM:
      length = "CM";
      break;
    default:
      length = m_length_unit;
  }

  m_event->use_units(momentum,length);
}

} // namespace Tauolapp
