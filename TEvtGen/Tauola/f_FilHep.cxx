#include "f_FilHep.h"
#include "Tauola.h"
#include "Log.h"

namespace Tauolapp
{

// NOTE: Not executed by release examples
float sgn( float a){
  return a/fabs(a);
}

// NOTE: Not executed by release examples
void filhepmc_print_details(int * n, int * status, int * pdg_id,
                            int * mother_first, int * mother_last, 
                            int * daughter_first, int * daughter_last, 
                            float p4[4], float * p_inv_mass, bool * photos_flag){

  Log::RedirectOutput(Log::Info());
  std::cout<<"*******************"<<std::endl;
  std::cout<<"Particle: "<<*n<<std::endl;
  std::cout<<"Status: "<<*status<<std::endl;
  std::cout<<"PDG ID: "<<*pdg_id<<std::endl;
  std::cout<<"Mothers: "<<*mother_first<<"-"<<*mother_last<<std::endl;
  std::cout<<"Daughters: "<<*daughter_first<<"-"<<*daughter_last<<std::endl;
  std::cout<<"4 momentum: "<<p4[0]<<", "<<p4[1]<<", "<<p4[2]<<", "<<p4[3]<<std::endl;
  std::cout<<"mass: "<<*p_inv_mass<<std::endl;
  std::cout<<"Photos Flag: "<<*photos_flag<<std::endl;
  std::cout<<"*******************"<<std::endl;
  Log::RevertOutput();
}


void filhep_(int * n, int * status, int * pdg_id,
             int * mother_first, int * mother_last, 
             int * daughter_first, int * daughter_last, 
             float p4[4], float * p_inv_mass, bool * photos_flag){

  /**  filhepmc_print_details(n, status, pdg_id,
                              mother_first, mother_last, 
                              daughter_first, daughter_last, 
                              p4, p_inv_mass, photos_flag);**/

  const int TAU_POSITION = 1;

  //Convert relative index's given by tauola to absolute ones:
  int abs_n=DecayList::getAbsoluteIndex(*n);

  //Create a new particle
  TauolaParticle * tau_mother = DecayList::getParticle(TAU_POSITION);
  TauolaParticle * new_particle;

  //  filhepmc_get_vertex(float v4[4]);  // vertex information add it to createNewParticle ...

  new_particle=tau_mother->createNewParticle(*pdg_id,*status,*p_inv_mass,
                                             p4[0],p4[1],p4[2],p4[3] );

  //boost along Z direction (Z defined as tau boost dir from tauola)
  if(Tauola::isUsingDecayOneBoost())
  {
     Tauola::decayOneBoost(tau_mother,new_particle);
  }
  else
  {
    if(tau_mother->getP(TauolaParticle::Z_AXIS)>0)
      new_particle->boostAlongZ(tau_mother->getP(),tau_mother->getE());
    else
      new_particle->boostAlongZ(-tau_mother->getP(),tau_mother->getE());
  }


  //Get rotation angles for transformation to lab frame.
  /**  double theta = tau_mother->getRotationAngle(TauolaParticle::Y_AXIS);
  tau_mother->rotate(TauolaParticle::Y_AXIS,theta);
  double phi = tau_mother->getRotationAngle(TauolaParticle::X_AXIS);
  tau_mother->rotate(TauolaParticle::Y_AXIS,-theta);

  //rotate coordinate system to lab frame.
  new_particle->rotate(TauolaParticle::X_AXIS,-phi);
  new_particle->rotate(TauolaParticle::Y_AXIS,-theta);**/

  //Add to list
  DecayList::updateList(new_particle, abs_n);

  //Get vector of mothers as TauolaParticles
  vector<TauolaParticle *> mothers;
  for(int i=*mother_first; i <= *mother_last && *mother_first!=0; i++){
    i=DecayList::getAbsoluteIndex(i,abs_n);
    mothers.push_back(DecayList::getParticle(i));
  }

  //Get vector of daughters as TauolaParticles
  vector<TauolaParticle *> daughters;
  for(int i=*daughter_first; i <= *daughter_last && *daughter_first!=0; i++){

    // NOTE: Not executed by release examples
    //       because daughter_first is always equal to 0
    i=DecayList::getAbsoluteIndex(i,abs_n);
    daughters.push_back(DecayList::getParticle(i));
  }

  //Add particle to event structure
  new_particle->setMothers(mothers);
  new_particle->setDaughters(daughters);

}

/** Simplified defintion. Only calculates mass (ams) from 4 momentum(p) */
void tralo4_(float * kto, float p[4], float q[4], float * ams){

  float tmp = p[3]*p[3] - p[1]*p[1] - p[2]*p[2] - p[0]*p[0];

  if (tmp!=0.0) tmp = tmp/sqrt(fabs(tmp));

  *ams = tmp;
}

} // namespace Tauolapp
