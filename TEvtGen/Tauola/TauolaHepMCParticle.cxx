#include "TauolaHepMCParticle.h"
#include "Log.h"

namespace Tauolapp
{

TauolaHepMCParticle::TauolaHepMCParticle(){
  m_particle = new HepMC::GenParticle();
}

TauolaHepMCParticle::~TauolaHepMCParticle(){
  
  //delete the mother and daughter pointers
  while(m_mothers.size()!=0){
    TauolaParticle * temp = m_mothers.back();
    m_mothers.pop_back();
    delete temp;
  }
  while(m_daughters.size()!=0){
    TauolaParticle * temp = m_daughters.back();
    m_daughters.pop_back();
    delete temp;
  }

  while(m_created_particles.size()!=0){
    TauolaHepMCParticle * temp = (TauolaHepMCParticle*) m_created_particles.back();
    m_created_particles.pop_back();
    if(temp->getHepMC()->barcode()==0) delete temp->getHepMC();
    delete temp;
  } 

}

// NOTE: Not executed by release examples
TauolaHepMCParticle::TauolaHepMCParticle(int pdg_id, int status, double mass){
  m_particle = new HepMC::GenParticle();
  m_particle->set_pdg_id(pdg_id);
  m_particle->set_status(status);
  m_particle->set_generated_mass(mass);
}

TauolaHepMCParticle::TauolaHepMCParticle(HepMC::GenParticle * particle){
  m_particle = particle;
}

HepMC::GenParticle * TauolaHepMCParticle::getHepMC(){
  return m_particle;
}

void TauolaHepMCParticle::undecay(){
  std::vector<TauolaParticle*> daughters = getDaughters();
  std::vector<TauolaParticle*>::iterator dIter = daughters.begin();

  for(; dIter != daughters.end(); dIter++)
    (*dIter)->undecay();

  if(m_particle->end_vertex())
  {
  while(m_particle->end_vertex()->particles_out_size())
  {
    HepMC::GenParticle *p = m_particle->end_vertex()->remove_particle(*(m_particle->end_vertex()->particles_out_const_begin()));
    delete p;
  }
  delete m_particle->end_vertex();
  }

  m_daughters.clear();
  m_particle->set_status(TauolaParticle::STABLE);

  for(unsigned int i=0;i<daughters.size();i++)
    delete daughters[i];
}

void TauolaHepMCParticle::setMothers(vector<TauolaParticle*> mothers){

  /******** Deal with mothers ***********/

  //If there are mothers
  if(mothers.size()>0){

    HepMC::GenParticle * part;
    part=dynamic_cast<TauolaHepMCParticle*>(mothers.at(0))->getHepMC();

    //Use end vertex of first mother as production vertex for particle
    HepMC::GenVertex * production_vertex = part->end_vertex();
    HepMC::GenVertex * orig_production_vertex = production_vertex;

    //If production_vertex does not exist - create it
    //If it's tau decay - set the time and position including the tau lifetime correction
    //otherwise - copy the time and position of decaying particle
    if(!production_vertex){
      production_vertex = new HepMC::GenVertex();
      HepMC::FourVector point = part->production_vertex()->position();
      production_vertex->set_position(point);
      part->parent_event()->add_vertex(production_vertex);
    }

    //Loop over all mothers to check that the end points to the right place
    vector<TauolaParticle*>::iterator mother_itr;
    for(mother_itr = mothers.begin(); mother_itr != mothers.end();
        mother_itr++){

      HepMC::GenParticle * moth;
      moth = dynamic_cast<TauolaHepMCParticle*>(*mother_itr)->getHepMC();

      if(moth->end_vertex()!=orig_production_vertex)
        Log::Fatal("Mother production_vertices point to difference places. Can not override. Please delete vertices first.",1);
      else
        production_vertex->add_particle_in(moth);

      //update status info
      if(moth->status()==TauolaParticle::STABLE)
        moth->set_status(TauolaParticle::DECAYED);
    }
    production_vertex->add_particle_out(m_particle);
  }
}

void TauolaHepMCParticle::setDaughters(vector<TauolaParticle*> daughters){

  if(!m_particle->parent_event())
    Log::Fatal("New particle needs the event set before it's daughters can be added",2);

  //If there are daughters
  if(daughters.size()>0){
    // NOTE: Not executed by release examples
    //       because daughters.size() is always 0

    //Use production vertex of first daughter as end vertex for particle
    HepMC::GenParticle * first_daughter;
    first_daughter = (dynamic_cast<TauolaHepMCParticle*>(daughters.at(0)))->getHepMC();

    HepMC::GenVertex * end_vertex;
    end_vertex=first_daughter->production_vertex();
    HepMC::GenVertex * orig_end_vertex = end_vertex;

    if(!end_vertex){ //if it does not exist create it
      end_vertex = new HepMC::GenVertex();
      m_particle->parent_event()->add_vertex(end_vertex);
    }

    //Loop over all daughters to check that the end points to the right place
    vector<TauolaParticle*>::iterator daughter_itr;
    for(daughter_itr = daughters.begin(); daughter_itr != daughters.end();
        daughter_itr++){

      HepMC::GenParticle * daug;
      daug = dynamic_cast<TauolaHepMCParticle*>(*daughter_itr)->getHepMC();


      if(daug->production_vertex()!=orig_end_vertex)
        Log::Fatal("Daughter production_vertices point to difference places. Can not override. Please delete vertices first.",3);
      else
        end_vertex->add_particle_out(daug);
    }
    end_vertex->add_particle_in(m_particle);
  }
}

std::vector<TauolaParticle*> TauolaHepMCParticle::getMothers(){

  if(m_mothers.size()==0&&m_particle->production_vertex()){
    HepMC::GenVertex::particles_in_const_iterator pcle_itr;
    pcle_itr=m_particle->production_vertex()->particles_in_const_begin();
    
    HepMC::GenVertex::particles_in_const_iterator pcle_itr_end;
    pcle_itr_end=m_particle->production_vertex()->particles_in_const_end();
    
    for(;pcle_itr != pcle_itr_end; pcle_itr++){
      m_mothers.push_back(new TauolaHepMCParticle(*pcle_itr));
    }
  }
  return m_mothers;
}

std::vector<TauolaParticle*> TauolaHepMCParticle::getDaughters(){
  
  if(m_daughters.size()==0&&m_particle->end_vertex()){
    HepMC::GenVertex::particles_out_const_iterator pcle_itr;
    pcle_itr=m_particle->end_vertex()->particles_out_const_begin();
    
    HepMC::GenVertex::particles_out_const_iterator pcle_itr_end;
    pcle_itr_end=m_particle->end_vertex()->particles_out_const_end();
    
    for(;pcle_itr != pcle_itr_end; pcle_itr++){
      m_daughters.push_back(new TauolaHepMCParticle(*pcle_itr));
    }
  }
  return m_daughters;
}

void TauolaHepMCParticle::checkMomentumConservation(){

  if(!m_particle->end_vertex()) return;
  
  // HepMC version of check_momentum_conservation
  // with added energy check

  double sumpx = 0, sumpy = 0, sumpz = 0, sume = 0;
  for( HepMC::GenVertex::particles_in_const_iterator part1 = m_particle->end_vertex()->particles_in_const_begin();
       part1 != m_particle->end_vertex()->particles_in_const_end(); part1++ ){

    sumpx += (*part1)->momentum().px();
    sumpy += (*part1)->momentum().py();
    sumpz += (*part1)->momentum().pz();
    sume  += (*part1)->momentum().e();
  }
  
  for( HepMC::GenVertex::particles_out_const_iterator part2 = m_particle->end_vertex()->particles_out_const_begin();
       part2 != m_particle->end_vertex()->particles_out_const_end(); part2++ ){

    sumpx -= (*part2)->momentum().px();
    sumpy -= (*part2)->momentum().py();
    sumpz -= (*part2)->momentum().pz();
    sume  -= (*part2)->momentum().e();
  }

  if( sqrt( sumpx*sumpx + sumpy*sumpy + sumpz*sumpz + sume*sume) > Tauola::momentum_conservation_threshold ) {
    Log::Warning()<<"Momentum not conserved in the vertex:"<<endl;
    Log::RedirectOutput(Log::Warning(false));
    m_particle->end_vertex()->print();
    Log::RevertOutput();
    return;
  }
  
  return;
}

// NOTE: Not executed by release examples
void TauolaHepMCParticle::setPdgID(int pdg_id){
  m_particle->set_pdg_id(pdg_id);
}

void TauolaHepMCParticle::setMass(double mass){
  m_particle->set_generated_mass(mass);
}

// NOTE: Not executed by release examples
void TauolaHepMCParticle::setStatus(int status){
  m_particle->set_status(status);
}

int TauolaHepMCParticle::getPdgID(){
  return m_particle->pdg_id();
}

int TauolaHepMCParticle::getStatus(){
  return m_particle->status();
}

int TauolaHepMCParticle::getBarcode(){
  return m_particle->barcode();
}

// Set (X,T) Position of tau decay trees
void TauolaHepMCParticle::decayEndgame(){

  double lifetime = Tauola::tau_lifetime * (-log( Tauola::randomDouble() ));
  HepMC::FourVector tau_momentum = m_particle->momentum();

  double mass     = sqrt(abs(  tau_momentum.e()*tau_momentum.e()
                             - tau_momentum.px()*tau_momentum.px()
                             - tau_momentum.py()*tau_momentum.py()
                             - tau_momentum.pz()*tau_momentum.pz()
                            ) );

  // Get previous position
  HepMC::FourVector previous_position = m_particle->production_vertex()->position();

  // Calculate new position
  HepMC::FourVector new_position(previous_position.x()+tau_momentum.px()/mass*lifetime,
                                 previous_position.y()+tau_momentum.py()/mass*lifetime,
                                 previous_position.z()+tau_momentum.pz()/mass*lifetime,
                                 previous_position.t()+tau_momentum.e() /mass*lifetime);

  // Set new position
  m_particle->end_vertex()->set_position(new_position);
  recursiveSetPosition(m_particle,new_position);
}

void TauolaHepMCParticle::recursiveSetPosition(HepMC::GenParticle *p, HepMC::FourVector pos){

  if(!p->end_vertex()) return;

  // Iterate over all outgoing particles
  for(HepMC::GenVertex::particles_out_const_iterator pp = p->end_vertex()->particles_out_const_begin();
      pp != p->end_vertex()->particles_out_const_end();
      ++pp){
    if( !(*pp)->end_vertex() ) continue;

    // Set position
    (*pp)->end_vertex()->set_position(pos);
    recursiveSetPosition(*pp,pos);
  }
}

TauolaHepMCParticle * TauolaHepMCParticle::createNewParticle(
                        int pdg_id, int status, double mass,
                        double px, double py, double pz, double e){

  TauolaHepMCParticle * new_particle = new TauolaHepMCParticle();
  new_particle->getHepMC()->set_pdg_id(pdg_id);
  new_particle->getHepMC()->set_status(status);
  new_particle->getHepMC()->set_generated_mass(mass);

  HepMC::FourVector momentum(px,py,pz,e);
  new_particle->getHepMC()->set_momentum(momentum);

  m_created_particles.push_back(new_particle);
  
  return new_particle;
}

void TauolaHepMCParticle::print(){
  m_particle->print();
}


/******** Getter and Setter methods: ***********************/

inline double TauolaHepMCParticle::getPx(){
  return m_particle->momentum().px();
}

inline double TauolaHepMCParticle::getPy(){
  return m_particle->momentum().py();
}

double TauolaHepMCParticle::getPz(){
  return m_particle->momentum().pz();
}

double TauolaHepMCParticle::getE(){
  return m_particle->momentum().e();
}

void TauolaHepMCParticle::setPx(double px){
  //make new momentum as something is wrong with
  //the HepMC momentum setters

  HepMC::FourVector momentum(m_particle->momentum());
  momentum.setPx(px);
  m_particle->set_momentum(momentum);
}

void TauolaHepMCParticle::setPy(double py){
  HepMC::FourVector momentum(m_particle->momentum());
  momentum.setPy(py);
  m_particle->set_momentum(momentum);
}


void TauolaHepMCParticle::setPz(double pz){
  HepMC::FourVector momentum(m_particle->momentum());
  momentum.setPz(pz);
  m_particle->set_momentum(momentum);
}

void TauolaHepMCParticle::setE(double e){
  HepMC::FourVector momentum(m_particle->momentum());
  momentum.setE(e);
  m_particle->set_momentum(momentum);
}

} // namespace Tauolapp
