#include "Tauola.h"
#include "TauolaParticlePair.h"
#include "Log.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>

namespace Tauolapp
{

/** constructor. Get the mothers, grandmothers and siblings of the tau */
TauolaParticlePair::TauolaParticlePair(std::vector<TauolaParticle*> &particle_list){
  TauolaParticle *particle=particle_list.back();
  particle_list.pop_back();
  setBornKinematics(0,0,-1.,0.); //set the default born variables

  // In case of decayOne() we need only the tau information
  if(Tauola::isUsingDecayOne())
  {
    m_mother = 0;
    m_mother_exists = false;
    m_production_particles.push_back(particle);
    m_final_particles.push_back(particle);
    initializeDensityMatrix();
    Log::AddDecay(0);
    return;
  }

  //See if there are any mothers
  std::vector<TauolaParticle *> temp_mothers = particle->findProductionMothers();

  // Case that there are no mothers or grandmothers
  if(temp_mothers.size()==0){
    // NOTE: Not executed by release examples
    //       However, such cases were present if tests used older Pythia8.1 or so.
    Log::Warning()<< "WARNING: Could not find taus mother or grandmothers. "
                  << "Ignoring spin effects" << std::endl;
    m_mother = 0;
    m_mother_exists = false;
    m_production_particles.push_back(particle);
    m_final_particles.push_back(particle->findLastSelf());
    initializeDensityMatrix();
    Log::AddDecay(1);
    return;
  }

  //fill the sibling pointers
  std::vector<TauolaParticle*> temp_daughters;
  temp_daughters = temp_mothers.at(0)->getDaughters();
  if(temp_daughters.size()==0)
    Log::Fatal("WARNING: Something wrong with event structure or there is a bug in the TAUOLA interface.",6);

  m_production_particles=temp_daughters;
  m_final_particles.push_back(particle->findLastSelf());

  //Second tau-like particle selection should match properties of the hard (exotic) process. At present, solution
  //match one of the possible diagrams contributing to SM process. Rare  case like tau+ tau+ nutau nutau will not be treted
  //accordingly to any diagram of SM
  for(signed int i=0; i < (int) m_production_particles.size(); i++)
  {
    //check if it has opposite PDGID sign
    if(m_final_particles.at(0)->getPdgID()*m_production_particles.at(i)->getPdgID()>=0) continue;

    if(m_production_particles.at(i)->getPdgID()==TauolaParticle::TAU_PLUS||
       m_production_particles.at(i)->getPdgID()==TauolaParticle::TAU_MINUS)
    {
       //tau+ or tau- - check if it's on the particle_list
       int j=-1;
       for(j=0;j<(int)particle_list.size();j++)
          if(m_production_particles.at(i)->getBarcode()==particle_list.at(j)->getBarcode()) break;
       if(j>=(int)particle_list.size()) continue;

       //exists on the list - add to m_final_particles and delete from particle_list
       m_final_particles.push_back(m_production_particles.at(i)->findLastSelf());
       particle_list.erase(particle_list.begin()+j);
    }
    else if(m_production_particles.at(i)->getPdgID()==TauolaParticle::TAU_NEUTRINO||
            m_production_particles.at(i)->getPdgID()==TauolaParticle::TAU_ANTINEUTRINO)
    {
       //neutrino - for now - just add to m_final_particles
       m_final_particles.push_back(m_production_particles.at(i)->findLastSelf());
    }
  }


  //fill the mother and grandmother pointers
  if(temp_mothers.size()==1){ //one mother
    m_mother_exists = true;
    m_mother = temp_mothers.at(0);
    m_grandmothers = m_mother->findProductionMothers();
    Log::AddDecay(3);
  }
  else{ //no mother, but grandparents exist
    m_mother_exists = false;
    m_grandmothers = temp_mothers;
    m_mother = makeTemporaryMother(m_production_particles);
    Log::AddDecay(2);
  } 

  initializeDensityMatrix();
  
  return;
}


/** The axis is defined by the boosting routine but our standard convention
    is:
    - Axis 3 is along the direction of the +ve tau, 
    - Axis 1 is perpendicular to the reaction plane (sign=??)
    - Axis 2 is defined through the vector product so the system
      is right handed (?? check). Axis 1,2 and 3 are parrellel for the 1st
      and second tau. */
void TauolaParticlePair::initializeDensityMatrix(){
      int incoming_pdg_id=0;
      int outgoing_pdg_id=0;
      double invariant_mass_squared=-5.0;
      double cosTheta=3.0;

  //initialize all elements of the density matrix to zero
  for(int x = 0; x < 4; x ++) {
    for(int y = 0; y < 4; y ++) 
      m_R[x][y] = 0;
  }

  m_R[0][0]=1;

  if(Tauola::isUsingDecayOne())
  {
    const double *pol = Tauola::getDecayOnePolarization();

    m_R[0][1]=pol[0];
    m_R[0][2]=pol[1];
    m_R[0][3]=pol[2];

    m_R[1][0]=pol[0];
    m_R[2][0]=pol[1];
    m_R[3][0]=pol[2];
  }

  if(!m_mother) return;
  // fill the matrix depending on mother


  // do scalar-pseudoscalar mixed higgs case separately since
  // switch doesn't allow non-constants in a case statement.
  // very annoying!

  if(m_mother->getPdgID()==Tauola::getHiggsScalarPseudoscalarPDG()){
    if(!Tauola::spin_correlation.HIGGS_H)  return;
    
    double phi = Tauola::getHiggsScalarPseudoscalarMixingAngle();
    double mass_ratio = Tauola::getTauMass()/m_mother->getMass();
    
    double beta = sqrt(1-4*mass_ratio*mass_ratio);
    double denominator = pow(cos(phi)*beta,2)+pow(sin(phi),2);
    
    m_R[0][0]=  1;
    m_R[1][1]=  (pow(cos(phi)*beta,2)-pow(sin(phi),2))/denominator;
    m_R[1][2]= 2*cos(phi)*sin(phi)*beta/denominator;
    m_R[2][1]= -m_R[1][2]; 
    m_R[2][2]=  m_R[1][1];
    m_R[3][3]= -1;
  }
  else {

    double pz = 0.0;

    switch(m_mother->getPdgID()){
      
    case TauolaParticle::Z0:
      if(!Tauola::spin_correlation.Z0)          break;
      // Here we calculate SVAR and COSTHE as well as IDE and IDF
      //   ANGULU(&IDE,&IDF,&SVAR,&COSTHE);
      // this is ++ configuration of Fig 5 from HEPPH0101311: 
      //pol_tau_minus[2]=-1; pol_tau_plus[2]=1;

      pz = getZPolarization(&incoming_pdg_id, &outgoing_pdg_id, &invariant_mass_squared, &cosTheta);
      m_R[0][0]=1;
      m_R[0][3]=2*pz-1;
      m_R[3][0]=2*pz-1;
      m_R[3][3]=1;
      // alternatively this may be overwritten if better solution exist
      recalculateRij(incoming_pdg_id, outgoing_pdg_id, invariant_mass_squared, cosTheta);
	  break;

    case TauolaParticle::GAMMA:
      if(!Tauola::spin_correlation.GAMMA)       break;
      // Here we calculate SVAR and COSTHE as well as IDE and IDF
      //   ANGULU(&IDE,&IDF,&SVAR,&COSTHE);
      // this is ++ configuration of Fig 5 from HEPPH0101311: 
      //pol_tau_minus[2]=-1; pol_tau_plus[2]=1;

      pz = getZPolarization(&incoming_pdg_id, &outgoing_pdg_id, &invariant_mass_squared, &cosTheta);
      m_R[0][0]=1;
      m_R[0][3]=2*pz-1;
      m_R[3][0]=2*pz-1;
      m_R[3][3]=1;
      // alternatively this may be overwritten if better solution exist
      recalculateRij(incoming_pdg_id, outgoing_pdg_id, invariant_mass_squared, cosTheta);
      break;

      //scalar higgs
    case TauolaParticle::HIGGS:
      if(!Tauola::spin_correlation.HIGGS)       break;
      m_R[0][0]=1;
      m_R[1][1]=1;
      m_R[2][2]=1;
      m_R[3][3]=-1;
      break;

      //pseudoscalar higgs case
    case TauolaParticle::HIGGS_A:
      if(!Tauola::spin_correlation.HIGGS_A)     break;
      m_R[0][0]=1;
      m_R[1][1]=-1;
      m_R[2][2]=-1;
      m_R[3][3]=-1;
      break;

    case TauolaParticle::HIGGS_PLUS:
      if(!Tauola::spin_correlation.HIGGS_PLUS)  break;
      m_R[0][0]=1;
      m_R[0][3]=-1;
      m_R[3][0]=-1;
      break;

    case TauolaParticle::HIGGS_MINUS:
      if(!Tauola::spin_correlation.HIGGS_MINUS) break;
      m_R[0][0]=1;
      m_R[0][3]=-1;
      m_R[3][0]=-1;
      break;

    case TauolaParticle::W_PLUS: 
      if(!Tauola::spin_correlation.W_PLUS)       break;
      m_R[0][0]=1;
      m_R[0][3]=1; //tau minus (tau minus is on the -ve Z axis)
      m_R[3][0]=1; //tau plus
      break;

    case TauolaParticle::W_MINUS:
      if(!Tauola::spin_correlation.W_MINUS)      break;
      m_R[0][0]=1;
      m_R[0][3]=1; //tau minus (tau minus is on the -ve Z axis)
      m_R[3][0]=1; //tau plus
      break;

    //ignore spin effects when mother is unknown
    default:
      m_R[0][0]=1;
      break;
    }
  }  

}

/**************************************************************/
void TauolaParticlePair::setBornKinematics(int incoming_pdg_id, int outgoing_pdg_id, double invariant_mass_squared,double cosTheta){
  Tauola::buf_incoming_pdg_id=incoming_pdg_id;
  Tauola::buf_outgoing_pdg_id=outgoing_pdg_id;
  Tauola::buf_invariant_mass_squared=invariant_mass_squared;
  Tauola::buf_cosTheta=cosTheta;
  //cout<<"(TauolaParticlePair::Just to be sure:) "<<buf_incoming_pdg_id<<" "<<buf_outgoing_pdg_id<<" "<<buf_invariant_mass_squared<<" "<<buf_cosTheta<<endl;
  //  m_R[0][0] to be added in next step;
}

/**************************************************************/
double TauolaParticlePair::getZPolarization(int *incoming_pdg_id, int *outgoing_pdg_id, double *invariant_mass_squared,double *cosTheta){

  //defaults
  *incoming_pdg_id = TauolaParticle::ELECTRON;
  *cosTheta = 0 ;
  *outgoing_pdg_id = TauolaParticle::TAU_PLUS;
  *invariant_mass_squared = pow(m_mother->getMass(),2);
  setBornKinematics(*incoming_pdg_id, *outgoing_pdg_id, *invariant_mass_squared, *cosTheta); // store for debugging

  //TRIVIAL CASE:
  //if we don't know the incoming beams then
  //return the average z polarisation
  if(m_grandmothers.size()<2){
      Log::Warning()<<"Not enough mothers of Z to "
                    <<"calculate cos(theta) between "
                    <<"incoming about outgoing beam"
                    << endl;
    return 1-plzap0_(incoming_pdg_id,outgoing_pdg_id, 
                     invariant_mass_squared, cosTheta);
  }

  if(!getTauPlus(m_production_particles)||!getTauMinus(m_production_particles)){
    Log::Error()<<"tau+ or tau- not found in Z decay"<< endl;
    return 0;
  }

  //NOW CHECK FOR THE DIFFICULT EVENTS:
  //case f1 + f2 + f3 -> Z -> tau tau
  //case f1 + f2 -> Z + f3, Z-> tau tau  or  f1 + f2 -> tau tau f3
  //case f1 + f2 -> Z -> tau tau gamma 
  if(m_grandmothers.size()>2 || 
     (m_grandmothers.at(0)->getDaughters().size()>1 && m_mother_exists==true) || 
     m_production_particles.size() > 2){

    //make a vector of the extra grandmother particles
    vector<TauolaParticle*> extra_grandmothers;
    for(int i=0; i<(int) m_grandmothers.size(); i++){
      //  temp_grandmothers.push_back(m_grandmothers.at(i));
      if(m_grandmothers.at(i)!=getGrandmotherPlus(m_grandmothers)&&
         m_grandmothers.at(i)!=getGrandmotherMinus(m_grandmothers))
         extra_grandmothers.push_back(m_grandmothers.at(i));
    }
    
    //make a vector of the tau siblings
    //and copy all the production particle vector.
    vector<TauolaParticle*> extra_tau_siblings;
    vector<TauolaParticle*> temp_production_particles;
    for(int i=0; i<(int) m_production_particles.size(); i++){
      if(m_production_particles.at(i)!=getTauPlus(m_production_particles)&&
         m_production_particles.at(i)!=getTauMinus(m_production_particles))
         extra_tau_siblings.push_back(m_production_particles.at(i));
    }

    //make a vector of the Z's sibling
    vector<TauolaParticle*> extra_Z_siblings;
    for(int i=0; m_mother_exists && i<(int) m_grandmothers.at(0)->getDaughters().size(); i++){
      if(m_grandmothers.at(0)->getDaughters().at(i)->getPdgID()!=TauolaParticle::Z0)
         extra_Z_siblings.push_back(m_grandmothers.at(0)->getDaughters().at(i));
    }

    //make temporary particles for the effect beams
    //and copy into a vector
    std::vector<TauolaParticle*> effective_taus;
    effective_taus.push_back(getTauPlus(m_production_particles)->clone());
    effective_taus.push_back(getTauMinus(m_production_particles)->clone());

    //copy grandmothers into the m_grandmothers vector since we want this to be perminante.
    TauolaParticle * g1 = getGrandmotherPlus(m_grandmothers)->clone();
    TauolaParticle * g2 = getGrandmotherMinus(m_grandmothers)->clone();
    m_grandmothers.clear();
    m_grandmothers.push_back(g1);
    m_grandmothers.push_back(g2);

    //loop over extra grandmothers
    for(int i=0; i<(int) extra_grandmothers.size(); i++)
      addToBeam(extra_grandmothers.at(i),&m_grandmothers,0); 

    //loop over siblings to the Z
    for(int i=0; i<(int) extra_Z_siblings.size(); i++)
      addToBeam(extra_Z_siblings.at(i),0,&m_grandmothers); 

    //loop over siblings of the taus
    for(int i=0; i<(int) extra_tau_siblings.size() ; i++){
      if(m_mother_exists)
        addToBeam(extra_tau_siblings.at(i),&effective_taus,0);
      else
        addToBeam(extra_tau_siblings.at(i),&effective_taus,&m_grandmothers);
    }
    //And we are finish making the effective income and outgoing beams

    TauolaParticle * temp_mother = makeTemporaryMother(effective_taus);
    *invariant_mass_squared = pow(temp_mother->getMass(),2);

   //now we are ready to do the boosting, calculate the angle
    //between incoming and outgoing, and get the polarisation of the Z

    double angle1,angle2,angle3;
    boostFromLabToTauPairFrame(&angle1, &angle2, &angle3, temp_mother,
                               m_grandmothers, effective_taus);

    /*double theta1 = -getGrandmotherPlus(m_grandmothers)->getRotationAngle(TauolaParticle::Y_AXIS);
    double theta2 = -(M_PI+getGrandmotherMinus(m_grandmothers)->getRotationAngle(TauolaParticle::Y_AXIS));
    *cosTheta = cos((theta1+theta2)/2.0); //just average the angles for now.*/

    TauolaParticle *tM=getTauPlus(effective_taus);
    TauolaParticle *gM=getGrandmotherPlus(m_grandmothers);

    double costheta1=(tM->getPx()*gM->getPx()+tM->getPy()*gM->getPy()+tM->getPz()*gM->getPz())/
                 sqrt(tM->getPx()*tM->getPx()+tM->getPy()*tM->getPy()+tM->getPz()*tM->getPz())/
                 sqrt(gM->getPx()*gM->getPx()+gM->getPy()*gM->getPy()+gM->getPz()*gM->getPz());
    tM=getTauMinus(effective_taus);
    gM=getGrandmotherMinus(m_grandmothers);

    double costheta2=(tM->getPx()*gM->getPx()+tM->getPy()*gM->getPy()+tM->getPz()*gM->getPz())/
                 sqrt(tM->getPx()*tM->getPx()+tM->getPy()*tM->getPy()+tM->getPz()*tM->getPz())/
                 sqrt(gM->getPx()*gM->getPx()+gM->getPy()*gM->getPy()+gM->getPz()*gM->getPz());
    double sintheta1 = sqrt(1-costheta1*costheta1);
    double sintheta2 = sqrt(1-costheta2*costheta2);
    double avg = (costheta1*sintheta2+costheta2*sintheta1)/(sintheta1+sintheta2);

    *cosTheta = avg;

    *incoming_pdg_id = getGrandmotherPlus(m_grandmothers)->getPdgID();

    if(abs(*incoming_pdg_id)==21 || abs(*incoming_pdg_id)==22)
    {
        *incoming_pdg_id = -getGrandmotherMinus(m_grandmothers)->getPdgID();
        //cout<<"INFO:\tgluon pdg id changed!"<<endl;
    }

    if( abs(*incoming_pdg_id)> 8  &&
        abs(*incoming_pdg_id)!=11 &&
        abs(*incoming_pdg_id)!=13 &&
        abs(*incoming_pdg_id)!=15 )
    {
      Log::Error()<<"Second class disaster: incoming_pdg_id = "<<*incoming_pdg_id<<endl;
      
      // Return avarage Z polarization
      *incoming_pdg_id = TauolaParticle::ELECTRON;
      *cosTheta = 0 ;
      *outgoing_pdg_id = TauolaParticle::TAU_PLUS;
      
      return 1-plzap0_(incoming_pdg_id,outgoing_pdg_id, 
                       invariant_mass_squared, cosTheta);
    }

    boostFromTauPairToLabFrame(angle1, angle2, angle3, temp_mother,
                               m_grandmothers, effective_taus);
    setBornKinematics(*incoming_pdg_id, *outgoing_pdg_id, *invariant_mass_squared, *cosTheta); // store for debugging    
    return 1-plzap0_(incoming_pdg_id,outgoing_pdg_id, invariant_mass_squared, cosTheta);

  }
  //REGULAR CASE:
  //f1 + f2 -> Z -> tau+ tau - or f1 f2 -> tau+ tau-
  //This includes Z -> tau tau, tau -> gamma tau?

  double angle1,angle2,angle3;
  boostFromLabToTauPairFrame(&angle1, &angle2, &angle3,
                             m_mother,m_grandmothers,m_production_particles);
 
  TauolaParticle *tM=getTauPlus(m_production_particles);
  TauolaParticle *gM=getGrandmotherPlus(m_grandmothers);
  double costheta1=(tM->getPx()*gM->getPx()+tM->getPy()*gM->getPy()+tM->getPz()*gM->getPz())/
               sqrt(tM->getPx()*tM->getPx()+tM->getPy()*tM->getPy()+tM->getPz()*tM->getPz())/
               sqrt(gM->getPx()*gM->getPx()+gM->getPy()*gM->getPy()+gM->getPz()*gM->getPz());

  tM=getTauMinus(m_production_particles);
  gM=getGrandmotherMinus(m_grandmothers);

  double costheta2=(tM->getPx()*gM->getPx()+tM->getPy()*gM->getPy()+tM->getPz()*gM->getPz())/
               sqrt(tM->getPx()*tM->getPx()+tM->getPy()*tM->getPy()+tM->getPz()*tM->getPz())/
               sqrt(gM->getPx()*gM->getPx()+gM->getPy()*gM->getPy()+gM->getPz()*gM->getPz());

  double sintheta1 = sqrt(1-costheta1*costheta1);
  double sintheta2 = sqrt(1-costheta2*costheta2);

  double avg       = (costheta1*sintheta2+costheta2*sintheta1)/(sintheta1+sintheta2);

  *cosTheta = avg;

  *incoming_pdg_id = getGrandmotherPlus(m_grandmothers)->getPdgID();

  boostFromTauPairToLabFrame(angle1, angle2, angle3,
                             m_mother,m_grandmothers,m_production_particles);

  setBornKinematics(*incoming_pdg_id, *outgoing_pdg_id, *invariant_mass_squared, *cosTheta); // store for debugging
  return 1-plzap0_(incoming_pdg_id,outgoing_pdg_id, invariant_mass_squared, cosTheta);
      // return 0.5 - (-0.12 + 0.12*cosTheta)/2;
}

/** WHERE WE CALCULATE THE EFFECTIVE BEAMS **/
/** This is where we decide which particle should be added into which beam,
    add it and change the flavour if necessary. candidates_same
    are on the same side of the vertex as the particle. This is needed
    for negative the particle 4-momentum. **/
void TauolaParticlePair::addToBeam(TauolaParticle * pcle,
                                   std::vector<TauolaParticle*> *candidates_same,
                                   std::vector<TauolaParticle*> *candidates_opp){


  double s=0.,o=-10.;
  TauolaParticle * s_beam = NULL;
  TauolaParticle * o_beam = NULL;

  if(candidates_same){
    double s0 =1.0/getVirtuality(pcle,candidates_same->at(0),false);
    double s1 = 1.0/getVirtuality(pcle,candidates_same->at(1),false);
    if(s0>s1){
      s=s0;
      s_beam = candidates_same->at(0);
    }
    else{
      s=s1;
      s_beam = candidates_same->at(1);
    }
  }
  if(candidates_opp){ 

    double o0 =1.0/getVirtuality(pcle,candidates_opp->at(0),true);
    double o1 = 1.0/getVirtuality(pcle,candidates_opp->at(1),true);
    if(o0>o1){
      o=o0;
      o_beam = candidates_opp->at(0);
    }
    else{
      o=o1;
      o_beam = candidates_opp->at(1);
    }
  }

  if(s>o)
  {
    s_beam->add(pcle);
    int pdg1 = pcle->getPdgID();
    int pdg2 = s_beam->getPdgID();
    if(abs(pdg2)==15) return;
    if((abs(pdg2)==21 || abs(pdg2)==22) && abs(pdg1)<17 && abs(pdg1)!=10 && abs(pdg1)!=9) s_beam->setPdgID( pdg1);
  }
  else
  {
    o_beam->subtract(pcle);
    int pdg1 = pcle->getPdgID();
    int pdg2 = o_beam->getPdgID();
    if((abs(pdg2)==21 || abs(pdg2)==22) && abs(pdg1)<17 && abs(pdg1)!=10 && abs(pdg1)!=9) o_beam->setPdgID(-pdg1);
  }

  //should we also do something with the flavours??
  
} 

double TauolaParticlePair::getVirtuality(TauolaParticle * p1, TauolaParticle*p2, bool flip){

  //if one particle in a gluon and the other a lepton
  if((p1->getPdgID()==TauolaParticle::GLUON&&abs(p2->getPdgID())>10&&abs(p2->getPdgID())<20) ||
     (p2->getPdgID()==TauolaParticle::GLUON&&abs(p1->getPdgID())>10&&abs(p1->getPdgID())<20))
    return -1;

  //add some others...

  //otherwise we calculate the virtuality:
  double kp = p1->getE()*p2->getE() - p1->getPx()*p2->getPx() 
    - p1->getPy()*p2->getPy() - p1->getPz()*p2->getPz();
  if(flip) // Note that we have 1/(p1+p2)^2  or 1/(p1-p2)^2 and p2^2=0
    kp = kp - p1->getMass()*p1->getMass();
 
  double q = 1;
  if(p1->getPdgID()==TauolaParticle::GAMMA){
    if(abs(p2->getPdgID())==TauolaParticle::UP || 
       abs(p2->getPdgID())==TauolaParticle::CHARM ||
       abs(p2->getPdgID())==TauolaParticle::TOP)
      q=2.0/3.0;
    if(abs(p2->getPdgID())==TauolaParticle::DOWN || 
       abs(p2->getPdgID())==TauolaParticle::STRANGE ||
       abs(p2->getPdgID())==TauolaParticle::BOTTOM)
      q=1.0/3.0;
  }
  if(p2->getPdgID()==TauolaParticle::GAMMA){
    if(abs(p1->getPdgID())==TauolaParticle::UP || 
       abs(p1->getPdgID())==TauolaParticle::CHARM ||
       abs(p1->getPdgID())==TauolaParticle::TOP)
      q=2.0/3.0;
    if(abs(p1->getPdgID())==TauolaParticle::DOWN || 
       abs(p1->getPdgID())==TauolaParticle::STRANGE ||
       abs(p1->getPdgID())==TauolaParticle::BOTTOM)
      q=1.0/3.0;
  }

  return fabs(2*kp)/q;
}

/***********************************************************
 **  TauolaParticlePair::decayTauPair().
 **  Generate tau decay
 ************************************************************/
void TauolaParticlePair::decayTauPair(){

  //initalize h vectors in case the partner is not a tau
  double h_tau_minus[4]={2,0,0,0}; //2 used to compensate for maximum weight
  double h_tau_plus[4]={2,0,0,0};  //in the case when there is only 1 tau
    
  TauolaParticle * tau_minus = getTauMinus(m_final_particles);
  TauolaParticle * tau_plus = getTauPlus(m_final_particles);

  //now calculate the spin weight
  for(double weight=0; weight < Tauola::randomDouble();){
    //call tauola decay and get polarimetric vectors
    if(tau_minus){
      Tauola::redefineTauMinusProperties(tau_minus);
      tau_minus->decay();
      h_tau_minus[0]=1;
      h_tau_minus[1]=tau_minus->getPolarimetricX();
      h_tau_minus[2]=tau_minus->getPolarimetricY();
      h_tau_minus[3]=tau_minus->getPolarimetricZ();
    }
    if(tau_plus){
      Tauola::redefineTauPlusProperties(tau_plus);
      tau_plus->decay();
      h_tau_plus[0]=1;
      h_tau_plus[1]=tau_plus->getPolarimetricX();
      h_tau_plus[2]=tau_plus->getPolarimetricY();
      h_tau_plus[3]=tau_plus->getPolarimetricZ();
    }
    //    cout<<" tau? tau? "<<  tau_plus->getPdgID()<<"  "<<  tau_minus->getPdgID()<<endl; 
    //    cout<<" przy uzyciu rrrRRRrrrRRR "<< m_R[0][0]<<" " <<m_R[3][3]<<" " << m_R[0][3]<<" " <<m_R[3][0] <<endl;
    //calculate the weight
    weight=0;
    for(int i =0; i<4; i++){
      for(int j=0; j<4; j++)
        weight+=m_R[i][j]*h_tau_minus[i]*h_tau_plus[j];
    }
    weight = weight/4.0;
  }
  double wthel[4]={0};
  wthel[0]=(h_tau_minus[0]+h_tau_minus[3])*(h_tau_plus[0]+h_tau_plus[3])*(m_R[0][0]+m_R[0][3]+m_R[3][0]+m_R[3][3]);
  wthel[1]=(h_tau_minus[0]+h_tau_minus[3])*(h_tau_plus[0]-h_tau_plus[3])*(m_R[0][0]-m_R[0][3]+m_R[3][0]-m_R[3][3]);
  wthel[2]=(h_tau_minus[0]-h_tau_minus[3])*(h_tau_plus[0]+h_tau_plus[3])*(m_R[0][0]+m_R[0][3]-m_R[3][0]-m_R[3][3]);
  wthel[3]=(h_tau_minus[0]-h_tau_minus[3])*(h_tau_plus[0]-h_tau_plus[3])*(m_R[0][0]-m_R[0][3]-m_R[3][0]+m_R[3][3]);

  double rn=Tauola::randomDouble();
  if (rn>(wthel[0]+wthel[1]+wthel[2])/(wthel[0]+wthel[1]+wthel[2]+wthel[3])) Tauola::setHelicities(-1,-1);
  else if (rn>(wthel[0]+wthel[1])    /(wthel[0]+wthel[1]+wthel[2]+wthel[3])) Tauola::setHelicities(-1,+1);
  else if (rn>(wthel[0])             /(wthel[0]+wthel[1]+wthel[2]+wthel[3])) Tauola::setHelicities( 1,-1);
  else                                                                       Tauola::setHelicities( 1, 1);


   

  //boost into the frame used to define the density matrices
  //and add the tau's new daughters.

  double angle1,angle2,angle3;
  TauolaParticle * mum = makeTemporaryMother(m_final_particles);

  if(!Tauola::isUsingDecayOneBoost())
       boostFromLabToTauPairFrame(&angle1, &angle2, &angle3,
                                  mum,m_grandmothers,m_final_particles);

  //add the accepted decay into the event record
  if(tau_plus)
    tau_plus->addDecayToEventRecord();
  if(tau_minus)
    tau_minus->addDecayToEventRecord();

  if(!Tauola::isUsingDecayOneBoost())
       boostFromTauPairToLabFrame(angle1,angle2,angle3,
                                  mum,m_grandmothers,m_final_particles);

  // Apply final decay modification, once taus are in lab frame. 
  // At present 28.02.11,  vertex position shift (in HepMC) is implemented.
  // Thanks Sho Iwamoto for feedback 
  if(tau_plus)
    tau_plus->decayEndgame();
  if(tau_minus)
    tau_minus->decayEndgame();

}

/***********************************************************
 **  Below are the methods used for boosting and rotating 
 **  into another reference frame.
 ************************************************************/

/** Step 1. (Transformation A). Any modification to this method also requires a modification
    to the inverse method boostFromTauPairFrameToLab (transformation A^-1).*/
void TauolaParticlePair::boostFromLabToTauPairFrame(double * rotation_angle1, 
                                                    double * rotation_angle2,
                                                    double * rotation_angle3,
                                                    TauolaParticle * mother,
                                                    vector<TauolaParticle *> grandmothers,
                                                    vector<TauolaParticle *> taus
                                                    ){ //in positive z axis (+1) //or negative z axis (-1)

  /** boost all gradmothers and daughters (taus, neutrinos, etc,) 
       to the mothers rest frame */

  for(int i=0; i< (int) grandmothers.size(); i++)
    grandmothers.at(i)->boostToRestFrame(mother);
  //If taus.size()==1 then taus[1] is the same as mum. Disaster to be avoided.
  if(taus.size()!=1)
    for(int i=0; i< (int) taus.size(); i++){
      taus.at(i)->boostToRestFrame(mother);
      // NOTE: Following line has no effect in release examples
      //       because taus don't have daughters at this step
      taus.at(i)->boostDaughtersToRestFrame(mother);
    }

  /** rotate all particles so taus are on the z axis */
  if(getTauPlus(taus)){ //if there's a tau+ use this to get the rotation angle
    *rotation_angle1 = getTauPlus(taus)->getRotationAngle(TauolaParticle::Y_AXIS);
    rotateSystem(grandmothers,taus,*rotation_angle1,TauolaParticle::Y_AXIS);
    *rotation_angle2 = getTauPlus(taus)->getRotationAngle(TauolaParticle::X_AXIS);
    rotateSystem(grandmothers,taus,*rotation_angle2,TauolaParticle::X_AXIS);
  }
  else{ //otherwise use the tau-
    *rotation_angle1 = M_PI+getTauMinus(taus)->getRotationAngle(TauolaParticle::Y_AXIS);
    rotateSystem(grandmothers,taus,*rotation_angle1,TauolaParticle::Y_AXIS);
    *rotation_angle2 = M_PI+getTauMinus(taus)->getRotationAngle(TauolaParticle::X_AXIS);
    rotateSystem(grandmothers,taus,*rotation_angle2,TauolaParticle::X_AXIS);
  }

  //now rotate incoming beams so they are aligned in the y-z plane
  if(grandmothers.size()<1){ //if there are no grandmothers
    *rotation_angle3 = 0;
    return;
  }

  //the first grandmother will have a positive y component
  if(getGrandmotherPlus(grandmothers))
    *rotation_angle3 = getGrandmotherPlus(grandmothers)->getRotationAngle(TauolaParticle::X_AXIS,TauolaParticle::Y_AXIS);
  else
    *rotation_angle3 = M_PI+getGrandmotherMinus(grandmothers)->getRotationAngle(TauolaParticle::X_AXIS,TauolaParticle::Y_AXIS);

  rotateSystem(grandmothers,taus,*rotation_angle3,TauolaParticle::X_AXIS,TauolaParticle::Y_AXIS);

}

/** Reverses boostFromLabtoMotherFrame. The three rotation angle must be provided.
    Any modification to this would require a modification to boostFromLabToTauPairFrame
    since this is the inverse transformation (Step 2: A^-1). */
void TauolaParticlePair::boostFromTauPairToLabFrame(double rotation_angle1, 
                                                    double rotation_angle2,
                                                    double rotation_angle3,
                                                    TauolaParticle * mother,
                                                    vector<TauolaParticle *> grandmothers,
                                                    vector<TauolaParticle *> taus){


  if(mother==0) //if there are no mothers
    return;

 
  //rotate grand mothers back away from the y-z plane
  rotateSystem(grandmothers,taus,-rotation_angle3, TauolaParticle::X_AXIS,TauolaParticle::Y_AXIS);

  //rotate all so that taus are no longer on the z axis
  rotateSystem(grandmothers,taus,-rotation_angle2,TauolaParticle::X_AXIS);
  rotateSystem(grandmothers,taus,-rotation_angle1,TauolaParticle::Y_AXIS);
 

  /*** boost grandmothers and daughters (taus) back into the lab frame */
  for(int i=0; i< (int) grandmothers.size(); i++)
    grandmothers.at(i)->boostFromRestFrame(mother);
  //If taus.size()==1 then taus[1] is the same as mum. Disaster to be avoided.
  if(taus.size()!=1)
    for(int i=0; i< (int) taus.size(); i++){
      taus.at(i)->boostFromRestFrame(mother);
      taus.at(i)->boostDaughtersFromRestFrame(mother);
    }
}

void TauolaParticlePair::rotateSystem(vector<TauolaParticle *> grandmothers,
                                      vector<TauolaParticle *> taus,
                                      double theta, int axis, int axis2){
  for(int i=0; i< (int) grandmothers.size(); i++)
    grandmothers.at(i)->rotate(axis,theta,axis2);
  for(int i=0; i< (int) taus.size(); i++){
    taus.at(i)->rotate(axis,theta,axis2);
    taus.at(i)->rotateDaughters(axis,theta,axis2);
  }
}



/***********************************************************
   Some extra useful methods
************************************************************/
TauolaParticle * TauolaParticlePair::makeTemporaryMother(vector<TauolaParticle *> taus){

    //calculate mothers 4-momentum based on the daughters.
    double e=0;
    double px=0;
    double py=0;
    double pz=0;
    
    bool tau_plus = false;
    bool tau_minus = false;
    bool tau_neutrino = false;
    bool tau_antineutrino = false;

    for(int d=0; d < (int) taus.size(); d++){
      TauolaParticle * daughter = taus.at(d);
      e+=daughter->getE();
      px+=daughter->getPx();
      py+=daughter->getPy();
      pz+=daughter->getPz();
      switch(daughter->getPdgID()){
         case TauolaParticle::TAU_PLUS:
           tau_plus=true;
           break;
         case TauolaParticle::TAU_MINUS:
           tau_minus=true;
           break;
         case TauolaParticle::TAU_NEUTRINO:
           tau_neutrino=true;
           break;
         case TauolaParticle::TAU_ANTINEUTRINO:
           tau_antineutrino=true;
           break;
      }
    }
    double m = e*e-px*px-py*py-pz*pz;
    if(m<0)
      m= -sqrt(-m);
    else
      m=sqrt(m);

    //look for mothers type:
    int pdg=0;

    //Assume Z
    if(tau_plus&&tau_minus)
      pdg=TauolaParticle::Z0 ;

    //Assume W+
    if(tau_plus&&tau_neutrino)
      pdg=TauolaParticle::W_PLUS;

    //Assume W-
    if(tau_minus&&tau_antineutrino)
      pdg=TauolaParticle::W_MINUS;

    // now make the mother
    return taus.at(0)->createNewParticle(pdg,2,m,px,py,pz,e);
}

// see if "particle" is one of the final taus in this tau pair 
// (based on particle barcode)
// NOTE: Not executed by release examples
bool TauolaParticlePair::contains(TauolaParticle * particle){

  for(int i=0; i< (int) m_final_particles.size(); i++){
    if(m_final_particles.at(i)->getBarcode()==particle->getBarcode())
      return true;
  } 
  return false;
}

TauolaParticle * TauolaParticlePair::getTauMinus(vector<TauolaParticle*> taus){
  for(int i=0; i< (int) taus.size(); i++){
    if(taus.at(i)->getPdgID()==TauolaParticle::TAU_MINUS) 
      return taus.at(i);
  }
  return 0;
}

TauolaParticle * TauolaParticlePair::getTauPlus(vector<TauolaParticle*> taus){
  for(int i=0; i< (int) taus.size(); i++){
    if(taus.at(i)->getPdgID()==TauolaParticle::TAU_PLUS)
      return taus.at(i);
  }
  return 0;
}

TauolaParticle * TauolaParticlePair::getGrandmotherPlus(vector<TauolaParticle*> grandmothers){
  //choose based on energy if there are more than one??
  double e = -1e30;
  int index = -1;
  for(int i=0; i< (int) grandmothers.size(); i++){
    if( (grandmothers.at(i)->getPdgID()<0 && grandmothers.at(i)->getPdgID()>-9) ||
       grandmothers.at(i)->getPdgID()== TauolaParticle::POSITRON||
       grandmothers.at(i)->getPdgID()== TauolaParticle::MUON_PLUS){
      if(e<grandmothers.at(i)->getE()){
        e=grandmothers.at(i)->getE();
        index=i;
      }
    }
  }
  if(index==-1)
  {
    for(int i=0; i< (int) grandmothers.size(); i++)
    {
      if(grandmothers.at(i)->getPdgID()==21 || grandmothers.at(i)->getPdgID()==22)
      {
        index=i;
        e=grandmothers.at(i)->getE();
        break;
      }
    }
  }
  if(index==-1) index=0;
  if(e==0)
  {
    grandmothers.at(index)->print();
    return 0;
  }
  else
    return grandmothers.at(index);
}

TauolaParticle * TauolaParticlePair::getGrandmotherMinus(vector<TauolaParticle*> grandmothers){
  //choose based on energy if there are more than one??
  double e = -1e30;
  int index = -1;
  for(int i=0; i< (int) grandmothers.size(); i++){
    if( (grandmothers.at(i)->getPdgID()>0 && grandmothers.at(i)->getPdgID()<9) ||
       grandmothers.at(i)->getPdgID()== TauolaParticle::ELECTRON||
       grandmothers.at(i)->getPdgID()== TauolaParticle::MUON_MINUS){
      if(e<grandmothers.at(i)->getE()){
        e=grandmothers.at(i)->getE();
        index=i;
      }
    }
  }
  if(index==-1)
  {
    for(int i=(int) grandmothers.size()-1; i>=0 ; i--)
    {
      if(grandmothers.at(i)->getPdgID()==21||grandmothers.at(i)->getPdgID()==22)
      {
        index=i;
        e=grandmothers.at(i)->getE();
        break;
      }
    }
  }
  if(index==-1) index=0;
  if(e==0)
    return 0;
  else
    return grandmothers.at(index);
}


   
void TauolaParticlePair::print(){

  Log::RedirectOutput(Log::Info());
  std::cout << "Daughters final:" << std::endl;
  for(int i=0; i< (int) m_final_particles.size(); i++)
    m_final_particles.at(i)->print();


  std::cout << "Daughters at production:" << std::endl;
  for(int i=0; i< (int) m_production_particles.size(); i++)
    m_production_particles.at(i)->print();


  std::cout << "Mother particle: " << std::endl;
  if(m_mother)
    m_mother->print();

  std::cout << "Grandmother particles: " << std::endl;
  for(int i=0; i< (int) m_grandmothers.size(); i++)
    m_grandmothers.at(i)->print();

  Log::RevertOutput();
}


void TauolaParticlePair::checkMomentumConservation(){

  for(int i=0; i< (int) m_final_particles.size(); i++)
    m_final_particles.at(i)->checkMomentumConservation();

  for(int i=0; i< (int) m_production_particles.size(); i++)
    m_production_particles.at(i)->checkMomentumConservation();

  if(m_mother)
    m_mother->checkMomentumConservation();

  for(int i=0; i< (int) m_grandmothers.size(); i++)
    m_grandmothers.at(i)->checkMomentumConservation();

}

void TauolaParticlePair::recalculateRij(int incoming_pdg_id, int outgoing_pdg_id, double invariant_mass_squared, double cosTheta){

  if (abs(outgoing_pdg_id)!=15) 
  {
    Log::Warning()<<"interface was not used for taus pdg id="<<outgoing_pdg_id<<endl;
    return;
  }

  double (*T) [Tauola::NCOS][4][4] = NULL;
  double (*Tw) [Tauola::NCOS]      = NULL;
  double (*Tw0)[Tauola::NCOS]      = NULL;
  double smin = 0.0;                       // Tauola::sminA;
  double smax = 0.0;                       // Tauola::smaxA;
  double step = 0.0;                       // (smaxl-sminl)/(Tauola::NS1-1);

  // Select table for appropriate incoming particles
  switch(abs(incoming_pdg_id)){
  case 11: 
    if(invariant_mass_squared<Tauola::smaxB && invariant_mass_squared>Tauola::sminB)
    { 
      T    = Tauola::table11B; 
      Tw   = Tauola::wtable11B; 
      Tw0  = Tauola::w0table11B; 
      smin = Tauola::sminB; 
      smax = Tauola::smaxB;
      step = (smax-smin)/(Tauola::NS2-1); 
    }
    else if (invariant_mass_squared<Tauola::smaxC && invariant_mass_squared>Tauola::sminC)
    { 
      T    = Tauola::table11C; 
      Tw   = Tauola::wtable11C; 
      Tw0  = Tauola::w0table11C; 
      smin = Tauola::sminC; 
      smax = Tauola::smaxC;
      step = (smax-smin)/(Tauola::NS3-1); 
    }
    else
    { 
      T    = Tauola::table11A; 
      Tw   = Tauola::wtable11A; 
      Tw0  = Tauola::w0table11A; 
      smin = Tauola::sminA; 
      smax = Tauola::smaxA;
      step = (smax-smin)/(Tauola::NS1-1); 
    }
    break;

  // NOTE: Use the same tables for 1, 3 and 5
  case  1:
  case  3:
  case  5: 
    if(invariant_mass_squared<Tauola::smaxB && invariant_mass_squared>Tauola::sminB)
    { 
      T    = Tauola::table1B; 
      Tw   = Tauola::wtable1B; 
      Tw0  = Tauola::w0table1B; 
      smin = Tauola::sminB; 
      smax = Tauola::smaxB;
      step = (smax-smin)/(Tauola::NS2-1); 
    }
    else if (invariant_mass_squared<Tauola::smaxC && invariant_mass_squared>Tauola::sminC)
    { 
      T    = Tauola::table1C; 
      Tw   = Tauola::wtable1C; 
      Tw0  = Tauola::w0table1C; 
      smin = Tauola::sminC; 
      smax = Tauola::smaxC;
      step = (smax-smin)/(Tauola::NS3-1); 
    }
    else
    { 
      T    = Tauola::table1A; 
      Tw   = Tauola::wtable1A; 
      Tw0  = Tauola::w0table1A; 
      smin = Tauola::sminA; 
      smax = Tauola::smaxA;
      step = (smax-smin)/(Tauola::NS1-1); 
    }
    break;

  // NOTE: Use the same tables for 2 and 4
  case  2:
  case  4: 
    if(invariant_mass_squared<Tauola::smaxB && invariant_mass_squared>Tauola::sminB)
    { 
      T    = Tauola::table2B; 
      Tw   = Tauola::wtable2B; 
      Tw0  = Tauola::w0table2B; 
      smin = Tauola::sminB; 
      smax = Tauola::smaxB;
      step = (smax-smin)/(Tauola::NS2-1); 
    }
    else if (invariant_mass_squared<Tauola::smaxC && invariant_mass_squared>Tauola::sminC)
    { 
      T    = Tauola::table2C; 
      Tw   = Tauola::wtable2C;
      Tw0  = Tauola::w0table2C;
      smin = Tauola::sminC; 
      smax = Tauola::smaxC;
      step = (smax-smin)/(Tauola::NS3-1); 
    }
    else
    { 
      T    = Tauola::table2A; 
      Tw   = Tauola::wtable2A;
      Tw0  = Tauola::w0table2A;
      smin = Tauola::sminA; 
      smax = Tauola::smaxA;
      step = (smax-smin)/(Tauola::NS1-1); 
    }
    break;

  default: 
     Log::Warning()<<"interface was not used for proper beams pdg id="<<incoming_pdg_id<<endl; 
     return;
  }

  // Check if we have found a table for matrix recalculation
  if (T[0][0][0][0]<0.5) return;

  // If mass is too small
  if (invariant_mass_squared <= exp(Tauola::sminA)){
    double mta   = 1.777; 
    double cosf  = cosTheta;
    double s     = invariant_mass_squared;
    double sinf  = sqrt(1-cosf*cosf);
    double MM    = sqrt(4e0*mta*mta/s);
    double xnorm = 1+cosf*cosf + MM*MM*sinf*sinf;

    // Zero matrix
    for (int i=0;i<4;i++)
      for (int j=0;j<4;j++)
        m_R[i][j]=0.0;

    m_R[0][0] = (1+cosf*cosf + MM*MM*sinf*sinf)/xnorm;
    m_R[1][1] = (-(1- MM*MM)*sinf*sinf)/xnorm;
    m_R[2][2] = ( (1+ MM*MM)*sinf*sinf)/xnorm;
    m_R[3][3] = (1+cosf*cosf - MM*MM*sinf*sinf)/xnorm;
    m_R[2][3] = (2*MM*sinf*cosf)/xnorm;
    m_R[3][2] = (2*MM*sinf*cosf)/xnorm;

    // Get weights
    double w  = 1.;
    double w0 = 1.;

    if     (Tauola::wtable2A[0][0]>0 ) w=Tauola::wtable2A[0][0];
    else if(Tauola::wtable1A[0][0]>0 ) w=Tauola::wtable1A[0][0];
    else if(Tauola::wtable11A[0][0]>0) w=Tauola::wtable11A[0][0];

    if     (Tauola::wtable2A[0][0]>0 ) w0=Tauola::w0table2A[0][0];
    else if(Tauola::wtable1A[0][0]>0 ) w0=Tauola::w0table1A[0][0];
    else if(Tauola::wtable11A[0][0]>0) w0=Tauola::w0table11A[0][0];

    //    Tauola::setEWwt(w/w0);
    Tauola::setEWwt(w,w0);
    return;
  } // if(mass too small)

  int    x       = 0;
  double buf     = smin;
  double remnant = 0.0;

  // Interpolate s
  if(T==Tauola::table11A || T==Tauola::table1A || T== Tauola::table2A)
  {
    while(log(invariant_mass_squared) > buf){
      x++;
      buf += (step);
    }
    remnant = (log(invariant_mass_squared)-(buf-step))/step;
  }
  else
  {
    while(invariant_mass_squared > buf){
      x++;
      buf += step;
    }
    remnant = (invariant_mass_squared-(buf-step))/step;
  }

  double cmin     = -1.;
  int    y        = 0;
  double remnantc = 0.0;

  // Interpolate cosTheta
  buf = cmin+1./Tauola::NCOS;
  while(cosTheta > buf){
    y++;
    buf+=2./Tauola::NCOS;
  }

  remnantc = (cosTheta-(buf-2./Tauola::NCOS))/(2./Tauola::NCOS);

  // Special cases at edges - extrapolation
  bool isUsingExtrapolation = false;
  if(y >= Tauola::NCOS) { isUsingExtrapolation = true; y = Tauola::NCOS-1; }
  if(y == 0)            { isUsingExtrapolation = true; y = 1;              }

  // Bilinear interpolation
  double b11,b21,b12,b22;

  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){

      if(isUsingExtrapolation){
        if(y == 1){
          b11 = 2*T[x-1][0][i][j] - T[x-1][1][i][j];
          b21 = 2*T[x][0][i][j]   - T[x][1][i][j];
          b12 = T[x-1][0][i][j];
          b22 = T[x][0][i][j];
        }
        else{
          b11 = T[x-1][y][i][j];
          b21 = T[x][y][i][j];
          b12 = 2*T[x-1][y][i][j] - T[x-1][y-1][i][j];
          b22 = 2*T[x][y][i][j]   - T[x][y-1][i][j];
        }
      } // if(isUsingExtrapolation)
      else{
        b11 = T[x-1][y-1][i][j];
        b21 = T[x][y-1][i][j];
        b12 = T[x-1][y][i][j];
        b22 = T[x][y][i][j];
      }
      m_R[i][j] = b11*(1-remnant)*(1-remnantc) + b21*remnant*(1-remnantc)
                + b12*(1-remnant)*remnantc + b22*remnant*remnantc;
    } // for(j)
  } // for(i)

  // Calculate electroweak weights
  double w,w0;

  if(isUsingExtrapolation){
    if(y == 1){
      b11 = 2*Tw[x-1][0] - Tw[x-1][1];
      b21 = 2*Tw[x][0]   - Tw[x][1];
      b12 = Tw[x-1][0];
      b22 = Tw[x][0];
    }
    else
    {
      b11 = Tw[x-1][y];
      b21 = Tw[x][y];
      b12 = 2*Tw[x-1][y] - Tw[x-1][y-1];
      b22 = 2*Tw[x][y]   - Tw[x][y-1];
    }
  } // if(isUsingExtrapolation)
  else
  {
    b11 = Tw[x-1][y-1];
    b21 = Tw[x][y-1];
    b12 = Tw[x-1][y];
    b22 = Tw[x][y];
  }

  w = b11*(1-remnant)*(1-remnantc) + b21*remnant*(1-remnantc)
    + b12*(1-remnant)*remnantc + b22*remnant*remnantc;

  if(isUsingExtrapolation){
    if(y == 1){
      b11 = 2*Tw0[x-1][0] - Tw0[x-1][1];
      b21 = 2*Tw0[x][0]   - Tw0[x][1];
      b12 = Tw0[x-1][0];
      b22 = Tw0[x][0];
    }
    else{
      b11 = Tw0[x-1][y];
      b21 = Tw0[x][y];
      b12 = 2*Tw0[x-1][y] - Tw0[x-1][y-1];
      b22 = 2*Tw0[x][y]   - Tw0[x][y-1];
    }
  } // if (isUsingExtrapolation)
  else{
    b11 = Tw0[x-1][y-1];
    b21 = Tw0[x][y-1];
    b12 = Tw0[x-1][y];
    b22 = Tw0[x][y];
  }

  w0 = b11*(1-remnant)*(1-remnantc) + b21*remnant*(1-remnantc)
     + b12*(1-remnant)*remnantc     + b22*remnant*remnantc;

  Tauola::setEWwt(w,w0);
}

} // namespace Tauolapp
