#include "TauolaHEPEVTEvent.h"

#include "Log.h"

namespace Tauolapp
{

TauolaHEPEVTEvent::~TauolaHEPEVTEvent()
{
  for(unsigned int i=0;i<particle_list.size();i++) delete particle_list[i];
}

TauolaHEPEVTEvent::TauolaHEPEVTEvent() {}

void TauolaHEPEVTEvent::addParticle(TauolaHEPEVTParticle *p)
{
  p->setEvent(this);

  p->setBarcode(particle_list.size());
  particle_list.push_back(p);
}

TauolaHEPEVTParticle *TauolaHEPEVTEvent::getParticle(int i)
{
  if( i<0 || i>=(int)particle_list.size() ) return NULL;
  return particle_list[i];
}

int TauolaHEPEVTEvent::getParticleCount()
{
  return particle_list.size();
}

// we have conflict in names, looks for -pdg_id also...
std::vector<TauolaParticle*> TauolaHEPEVTEvent::findParticles(int pdg_id){

  std::vector<TauolaParticle*> list;

  // Loop over all particles in the event looking
  // for tau (or other) particle with specified pdg_id and -pdg_id
  for(unsigned int i=0; i<particle_list.size(); i++)
  {
    if( abs(particle_list[i]->getPdgID() ) == pdg_id)
      list.push_back(particle_list[i]);
  }

  return list;
}

// we have conflict in names, should be findStableTaus or have another argument. 
std::vector<TauolaParticle*> TauolaHEPEVTEvent::findStableParticles(int pdg_id){

  std::vector<TauolaParticle*> tau_list = findParticles(pdg_id);
  std::vector<TauolaParticle*> stable_tau_list;

  for(int i=0; i<(int) tau_list.size(); i++){

    if(!tau_list.at(i)->hasDaughters())
      stable_tau_list.push_back(tau_list[i]);
    else
    {
      std::vector<TauolaParticle*> t = tau_list[i]->getDaughters();
      //Ignore taus that we won't be decaying anyway
      if(t.size()==1) continue;
      if(t.size()==2 && (abs(t[0]->getPdgID())==15 || abs(t[1]->getPdgID())==15) ) continue;
      Log::Warning()<<"Particle with pdg code "<<tau_list.at(i)->getPdgID()
                    <<" already has daughters" <<endl;
    }
  }

  return stable_tau_list;

}

void TauolaHEPEVTEvent::print()
{
  printf("TauolaHEPEVTEvent\n-----------------\n");
  for(unsigned int i=0;i<particle_list.size();i++) particle_list[i]->print();
}

void TauolaHEPEVTEvent::clear()
{
  for(unsigned int i=0;i<particle_list.size();i++) delete particle_list[i];
  particle_list.clear();
}

#ifdef USE_HEPEVT_INTERFACE

void TauolaHEPEVTEvent::read_event_from_HEPEVT(TauolaHEPEVTEvent *evt)
{
  if(evt==NULL) return;
  
  for(int i=0; i<hepevt_.nhep; i++)
  {
    TauolaHEPEVTParticle *p = new TauolaHEPEVTParticle
    (
      hepevt_.idhep [i],
      hepevt_.isthep[i],
      hepevt_.phep  [i][0],
      hepevt_.phep  [i][1],
      hepevt_.phep  [i][2],
      hepevt_.phep  [i][3],
      hepevt_.phep  [i][4],
      hepevt_.jmohep[i][0]-1,
      hepevt_.jmohep[i][1]-1,
      hepevt_.jdahep[i][0]-1,
      hepevt_.jdahep[i][1]-1
    );
    evt->addParticle(p);
  }
}

void TauolaHEPEVTEvent::write_event_to_HEPEVT(TauolaHEPEVTEvent *evt)
{
  if(evt==NULL) return;
  
  hepevt_.nhep = evt->getParticleCount();
  
  for(int i=0; i<hepevt_.nhep; i++)
  {
    TauolaHEPEVTParticle *p = evt->getParticle(i);
    
    hepevt_.idhep [i]   =p->getPdgID();
    hepevt_.isthep[i]   =p->getStatus();
    hepevt_.phep  [i][0]=p->getPx();
    hepevt_.phep  [i][1]=p->getPy();
    hepevt_.phep  [i][2]=p->getPz();
    hepevt_.phep  [i][3]=p->getE();
    hepevt_.phep  [i][4]=p->getMass();
    hepevt_.jmohep[i][0]=p->getFirstMotherIndex()  +1;
    hepevt_.jmohep[i][1]=p->getSecondMotherIndex() +1;
    hepevt_.jdahep[i][0]=p->getDaughterRangeStart()+1;
    hepevt_.jdahep[i][1]=p->getDaughterRangeEnd()  +1;
    hepevt_.vhep  [i][0]=0.0;
    hepevt_.vhep  [i][1]=0.0;
    hepevt_.vhep  [i][2]=0.0;
    hepevt_.vhep  [i][3]=0.0;
  }
}

#endif

} // namespace Tauolapp
