#include "DecayList.h"
#include "Log.h"

using namespace std;

namespace Tauolapp
{

vector<TauolaParticle*> DecayList::m_particle_list;

int DecayList::getAbsoluteIndex(int index){
  return getAbsoluteIndex(index, m_particle_list.size()+1);
}

int DecayList::getAbsoluteIndex(int index, 
                                int neg_index_relative_to){
  int absIndex;

  if(index > 0) //absolute position
    absIndex = index;
  else //relative to fixed
    absIndex = index + neg_index_relative_to;
  //Some error checking
  if(absIndex < 1 || absIndex > (int)m_particle_list.size()+1){
    Log::Error()<<"Index outside range: "<< absIndex << ". Range: 1 to "
                << m_particle_list.size()+1 << endl;
    Log::Fatal(4);
  }
  //  cout << "Final call in getAbsoluteIndex().. "<< absIndex << endl;
  return absIndex; //account for vectors starting at index 0
}

// NOTE: Not executed by release examples
int DecayList::getAbsoluteIndex(TauolaParticle * particle){
  for(int i=0; i < (int) m_particle_list.size(); i++){
    if(m_particle_list.at(i)==particle)
      return i+1;
  }
  Log::Warning()<<"Could not find particle in particle_list" << endl;
  return 0;
}

TauolaParticle * DecayList::getParticle(int index){
  return m_particle_list.at(index-1);
}

void DecayList::updateList(TauolaParticle * new_particle,
                           int index){
  
  if(index > (int) m_particle_list.size())
    //Add new particle to end
    addToEnd(new_particle);
  else{ 
    // NOTE: Not executed by release examples

    TauolaParticle * old_particle = getParticle(index);
    //Add new particle
    m_particle_list.at(index - 1) = new_particle;

    //Remove old particle at same index in event record
    /**    if(old_particle->production_vertex())
      old_particle->production_vertex()->remove_particle(old_particle);
    if(old_particle->end_vertex())
      old_particle->end_vertex()->remove_particle(old_particle);
      delete old_particle;**/
    delete old_particle;

  }
}

void DecayList::addToEnd(TauolaParticle * new_particle){
  m_particle_list.push_back(new_particle);
}

void DecayList::print(){
  for(int index=0; index < (int) m_particle_list.size(); index++){
    Log::Info()<< "Index: "<< index+1<<" Object: "<< m_particle_list.at(index)<<endl;
  }
}
  
void DecayList::clear(){
  m_particle_list.clear();
}

} // namespace Tauolapp
