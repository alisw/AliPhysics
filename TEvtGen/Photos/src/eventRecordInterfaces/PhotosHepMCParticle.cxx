#include "HepMC/GenEvent.h"
#include "PhotosHepMCParticle.h"
#include "Log.h"
#include "Photos.h"

namespace Photospp
{

PhotosHepMCParticle::PhotosHepMCParticle(){
  m_particle = new HepMC::GenParticle();
}

PhotosHepMCParticle::PhotosHepMCParticle(int pdg_id, int status, double mass){
  m_particle = new HepMC::GenParticle();
  m_particle->set_pdg_id(pdg_id);
  m_particle->set_status(status);
  m_particle->set_generated_mass(mass);
}

PhotosHepMCParticle::PhotosHepMCParticle(HepMC::GenParticle * particle){
  m_particle = particle;
}

PhotosHepMCParticle::~PhotosHepMCParticle(){
  clear(m_mothers);
  clear(m_daughters);
  //  clear(m_created_particles);
}
 

//delete the TauolaHepMCParticle objects
void PhotosHepMCParticle::clear(std::vector<PhotosParticle*> v){
  while(v.size()!=0){
    PhotosParticle * temp = v.back();
    v.pop_back();
    delete temp;
  }
}

HepMC::GenParticle * PhotosHepMCParticle::getHepMC(){
  return m_particle;
}

void PhotosHepMCParticle::setMothers(vector<PhotosParticle*> mothers){

  /******** Deal with mothers ***********/

  clear(m_mothers);

  //If there are mothers
  if(mothers.size()>0){

    HepMC::GenParticle * part;
    part=dynamic_cast<PhotosHepMCParticle*>(mothers.at(0))->getHepMC();

    //Use end vertex of first mother as production vertex for particle
    HepMC::GenVertex * production_vertex = part->end_vertex();
    HepMC::GenVertex * orig_production_vertex = production_vertex;

    if(!production_vertex){ //if it does not exist create it
      production_vertex = new HepMC::GenVertex();
      part->parent_event()->add_vertex(production_vertex);
    }

    //Loop over all mothers to check that the end points to the right place
    vector<PhotosParticle*>::iterator mother_itr;
    for(mother_itr = mothers.begin(); mother_itr != mothers.end();
        mother_itr++){

      HepMC::GenParticle * moth;
      moth = dynamic_cast<PhotosHepMCParticle*>(*mother_itr)->getHepMC();

      if(moth->end_vertex()!=orig_production_vertex)
        Log::Fatal("PhotosHepMCParticle::setMothers(): Mother production_vertices point to difference places. Can not override. Please delete vertices first.",1);
      else
        production_vertex->add_particle_in(moth);

      //update status info
      if(moth->status()==PhotosParticle::STABLE)
        moth->set_status(PhotosParticle::DECAYED);
    }
    production_vertex->add_particle_out(m_particle);
  }
}



void PhotosHepMCParticle::addDaughter(PhotosParticle* daughter){

  //add to this classes internal list as well.
  m_daughters.push_back(daughter);

  //this assumes there is already an end vertex for the particle

  if(!m_particle->end_vertex())
    Log::Fatal("PhotosHepMCParticle::addDaughter(): This method assumes an end_vertex exists. Maybe you really want to use setDaughters.",2);

  HepMC::GenParticle * daugh = (dynamic_cast<PhotosHepMCParticle*>(daughter))->getHepMC();
  m_particle->end_vertex()->add_particle_out(daugh);
  
}

void PhotosHepMCParticle::setDaughters(vector<PhotosParticle*> daughters){

  if(!m_particle->parent_event())
    Log::Fatal("PhotosHepMCParticle::setDaughters(): New particle needs the event set before it's daughters can be added",3);

  clear(m_daughters);

  //If there are daughters
  if(daughters.size()>0){

    //Use production vertex of first daughter as end vertex for particle
    HepMC::GenParticle * first_daughter;
    first_daughter = (dynamic_cast<PhotosHepMCParticle*>(daughters.at(0)))->getHepMC();

    HepMC::GenVertex * end_vertex;
    end_vertex=first_daughter->production_vertex();
    HepMC::GenVertex * orig_end_vertex = end_vertex;

    if(!end_vertex){ //if it does not exist create it
      end_vertex = new HepMC::GenVertex();
      m_particle->parent_event()->add_vertex(end_vertex);
    }

    //Loop over all daughters to check that the end points to the right place
    vector<PhotosParticle*>::iterator daughter_itr;
    for(daughter_itr = daughters.begin(); daughter_itr != daughters.end();
        daughter_itr++){

      HepMC::GenParticle * daug;
      daug = dynamic_cast<PhotosHepMCParticle*>(*daughter_itr)->getHepMC();


      if(daug->production_vertex()!=orig_end_vertex)
        Log::Fatal("PhotosHepMCParticle::setDaughters(): Daughter production_vertices point to difference places. Can not override. Please delete vertices first.",4);
      else
        end_vertex->add_particle_out(daug);
    }
    end_vertex->add_particle_in(m_particle);
  }

}

std::vector<PhotosParticle*> PhotosHepMCParticle::getMothers(){

  if(m_mothers.size()==0&&m_particle->production_vertex()){

    HepMC::GenVertex::particles_in_const_iterator pcle_itr;
    pcle_itr=m_particle->production_vertex()->particles_in_const_begin();
    
    HepMC::GenVertex::particles_in_const_iterator pcle_itr_end;
    pcle_itr_end=m_particle->production_vertex()->particles_in_const_end();
    
    for(;pcle_itr != pcle_itr_end; pcle_itr++){
      m_mothers.push_back(new PhotosHepMCParticle(*pcle_itr));
    }
  }
  return m_mothers;
}

std::vector<PhotosParticle*> PhotosHepMCParticle::getDaughters(){

  if(m_daughters.size()==0&&m_particle->end_vertex()){
    HepMC::GenVertex::particles_out_const_iterator pcle_itr;
    pcle_itr=m_particle->end_vertex()->particles_out_const_begin();
    
    HepMC::GenVertex::particles_out_const_iterator pcle_itr_end;
    pcle_itr_end=m_particle->end_vertex()->particles_out_const_end();
    
    for(;pcle_itr != pcle_itr_end; pcle_itr++){

      // ommit particles if their status code is ignored by Photos
      if( Photos::isStatusCodeIgnored( (*pcle_itr)->status() ) ) continue;

      m_daughters.push_back(new PhotosHepMCParticle(*pcle_itr));
    }
  }
  return m_daughters;

}

std::vector<PhotosParticle*> PhotosHepMCParticle::getAllDecayProducts(){

  m_decay_products.clear();

  if(!hasDaughters()) // if this particle has no daughters
    return m_decay_products;

  std::vector<PhotosParticle*> daughters = getDaughters();
  
  // copy daughters to list of all decay products
  m_decay_products.insert(m_decay_products.end(),daughters.begin(),daughters.end());
  
  // Now, get all daughters recursively, without duplicates.
  // That is, for each daughter:
  // 1)  get list of her daughters
  // 2)  for each particle on this list:
  //  a) check if it is already on the list
  //  b) if it's not, add her to the end of the list
  for(unsigned int i=0;i<m_decay_products.size();i++)
  {
    std::vector<PhotosParticle*> daughters2 = m_decay_products[i]->getDaughters();

    if(!m_decay_products[i]->hasDaughters()) continue;
    for(unsigned int j=0;j<daughters2.size();j++)
    {
      bool add=true;
      for(unsigned int k=0;k<m_decay_products.size();k++)
        if( daughters2[j]->getBarcode() == m_decay_products[k]->getBarcode() )
        {
          add=false;
          break;
        }

      if(add) m_decay_products.push_back(daughters2[j]);
    }
  }

  return m_decay_products;
}

bool PhotosHepMCParticle::checkMomentumConservation(){

  if(!m_particle->end_vertex()) return true;
  
  // HepMC version of check_momentum_conservation
  // Ommitting history entries (status == 3)

  double sumpx = 0, sumpy = 0, sumpz = 0, sume = 0;
  for( HepMC::GenVertex::particles_in_const_iterator part1 = m_particle->end_vertex()->particles_in_const_begin();
       part1 != m_particle->end_vertex()->particles_in_const_end(); part1++ ){

    if( Photos::isStatusCodeIgnored((*part1)->status()) ) continue;

    sumpx += (*part1)->momentum().px();
    sumpy += (*part1)->momentum().py();
    sumpz += (*part1)->momentum().pz();
    sume  += (*part1)->momentum().e();
  }
  
  for( HepMC::GenVertex::particles_out_const_iterator part2 = m_particle->end_vertex()->particles_out_const_begin();
       part2 != m_particle->end_vertex()->particles_out_const_end(); part2++ ){

    if( Photos::isStatusCodeIgnored((*part2)->status()) ) continue;

    sumpx -= (*part2)->momentum().px();
    sumpy -= (*part2)->momentum().py();
    sumpz -= (*part2)->momentum().pz();
    sume  -= (*part2)->momentum().e();
  }

  if( sqrt( sumpx*sumpx + sumpy*sumpy + sumpz*sumpz + sume*sume) > Photos::momentum_conservation_threshold ) {
    Log::Warning()<<"Momentum not conserved in the vertex:"<<endl;
    Log::RedirectOutput(Log::Warning(false));
    m_particle->end_vertex()->print();
    Log::RevertOutput();
    return false;
  }
  
  return true;
}

void PhotosHepMCParticle::setPdgID(int pdg_id){
  m_particle->set_pdg_id(pdg_id);
}

void PhotosHepMCParticle::setMass(double mass){
  m_particle->set_generated_mass(mass);
}

void PhotosHepMCParticle::setStatus(int status){
  m_particle->set_status(status);
}

int PhotosHepMCParticle::getPdgID(){
  return m_particle->pdg_id();
}

int PhotosHepMCParticle::getStatus(){
  return m_particle->status();
}

int PhotosHepMCParticle::getBarcode(){
  return m_particle->barcode();
}


PhotosHepMCParticle * PhotosHepMCParticle::createNewParticle(
                        int pdg_id, int status, double mass,
                        double px, double py, double pz, double e){

  PhotosHepMCParticle * new_particle = new PhotosHepMCParticle();
  new_particle->getHepMC()->set_pdg_id(pdg_id);
  new_particle->getHepMC()->set_status(status);
  new_particle->getHepMC()->set_generated_mass(mass);

  HepMC::FourVector momentum(px,py,pz,e);
  new_particle->getHepMC()->set_momentum(momentum);

  m_created_particles.push_back(new_particle);
  return new_particle;
}

void PhotosHepMCParticle::createHistoryEntry(){

  if(!m_particle->production_vertex())
  {
    Log::Warning()<<"PhotosHepMCParticle::createHistoryEntry(): particle without production vertex."<<endl;
    return;
  }
  
  HepMC::GenParticle *part = new HepMC::GenParticle(*m_particle);
  part->set_status(Photos::historyEntriesStatus);
  m_particle->production_vertex()->add_particle_out(part);
}

void PhotosHepMCParticle::createSelfDecayVertex(PhotosParticle *out)
{
  if(m_particle->end_vertex())
  {
    Log::Error()<<"PhotosHepMCParticle::createSelfDecayVertex: particle already has end vertex!"<<endl;
    return;
  }

  if(getHepMC()->parent_event()==NULL)
  {
    Log::Error()<<"PhotosHepMCParticle::createSelfDecayVertex: particle not in the HepMC event!"<<endl;
    return;
  }

  // Add new vertex and new particle to HepMC
  HepMC::GenParticle *outgoing = new HepMC::GenParticle( *(dynamic_cast<PhotosHepMCParticle*>(out)->m_particle) );
  HepMC::GenVertex   *v        = new HepMC::GenVertex();
  
  // Copy vertex position from parent vertex
  v->set_position( m_particle->production_vertex()->position() );
  
  v->add_particle_in (m_particle);
  v->add_particle_out(outgoing);
  
  getHepMC()->parent_event()->add_vertex(v);
  
  // If this particle was stable, set its status to 2
  if(getStatus()==1) setStatus(2);
}

void PhotosHepMCParticle::print(){
  m_particle->print();
}


/******** Getter and Setter methods: ***********************/

inline double PhotosHepMCParticle::getPx(){
  return m_particle->momentum().px();
}

inline double PhotosHepMCParticle::getPy(){
  return m_particle->momentum().py();
}

double PhotosHepMCParticle::getPz(){
  return m_particle->momentum().pz();
}

double PhotosHepMCParticle::getE(){
  return m_particle->momentum().e();
}

void PhotosHepMCParticle::setPx(double px){
  //make new momentum as something is wrong with
  //the HepMC momentum setters

  HepMC::FourVector momentum(m_particle->momentum());
  momentum.setPx(px);
  m_particle->set_momentum(momentum);
}

void PhotosHepMCParticle::setPy(double py){
  HepMC::FourVector momentum(m_particle->momentum());
  momentum.setPy(py);
  m_particle->set_momentum(momentum);
}


void PhotosHepMCParticle::setPz(double pz){
  HepMC::FourVector momentum(m_particle->momentum());
  momentum.setPz(pz);
  m_particle->set_momentum(momentum);
}

void PhotosHepMCParticle::setE(double e){
  HepMC::FourVector momentum(m_particle->momentum());
  momentum.setE(e);
  m_particle->set_momentum(momentum);
}

double PhotosHepMCParticle::getMass()
{
	return m_particle->generated_mass();
}

} // namespace Photospp
