#include "TauolaParticle.h"
#include "Log.h"
#include <stdexcept> 

namespace Tauolapp
{

double TauolaParticle::getPolarimetricX(){
  return m_pol_x;
}

double TauolaParticle::getPolarimetricY(){
  return m_pol_y;
}

double TauolaParticle::getPolarimetricZ(){
  return m_pol_z;
}


TauolaParticle * TauolaParticle::clone(){
  
  return createNewParticle(getPdgID(),getStatus(),getMass(),
                           getPx(),getPy(),getPz(),getE());

}

// NOTE: Not executed by release examples
double TauolaParticle::getAngle(TauolaParticle * other_particle){

  //use the dot product
  double x1 = getPx();
  double y1 = getPy();
  double z1 = getPz();
  double x2 = other_particle->getPx();
  double y2 = other_particle->getPy();
  double z2 = other_particle->getPz();

  return acos( (x1*x2+y1*y2+z1*z2) / sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2)) );

}

// NOTE: Not executed by release examples
void TauolaParticle::add(TauolaParticle * other_particle){

  setPx(getPx() + other_particle->getPx());
  setPy(getPy() + other_particle->getPy());
  setPz(getPz() + other_particle->getPz());
  setE(getE() + other_particle->getE());
  setMass( sqrt( getE()*getE()-getPx()*getPx()-getPy()*getPy()-getPz()*getPz() ));
}

void TauolaParticle::subtract(TauolaParticle * other_particle){

  setPx(getPx() - other_particle->getPx());
  setPy(getPy() - other_particle->getPy());
  setPz(getPz() - other_particle->getPz());
  setE(getE() - other_particle->getE());
  setMass( sqrt( getE()*getE()-getPx()*getPx()-getPy()*getPy()-getPz()*getPz() ));
}

int TauolaParticle::getSign(){
  if(getPdgID()== Tauola::getDecayingParticle())
    return SAME_SIGN;
  else if(getPdgID()== -1 * Tauola::getDecayingParticle())
    return OPPOSITE_SIGN;
  else
    return NA_SIGN;
}

bool TauolaParticle::hasDaughters(){
  if(getDaughters().size()==0)
    return 0;
  else
    return 1;
}

TauolaParticle * TauolaParticle::findLastSelf(){
  vector<TauolaParticle*> daughters = getDaughters();
  vector<TauolaParticle*>::iterator pcl_itr = daughters.begin();
  
  //get all daughters and look for stable with same pgd id
  for(;pcl_itr != daughters.end();pcl_itr++){
    if((*pcl_itr)->getPdgID()==this->getPdgID())
      return (*pcl_itr)->findLastSelf();
  }

  return this;
}

std::vector<TauolaParticle*> TauolaParticle::findProductionMothers(){
  vector<TauolaParticle*> mothers = getMothers();
  vector<TauolaParticle*>::iterator pcl_itr = mothers.begin();
  
    //get all mothers and check none have pdg id of this one
    for(;pcl_itr != mothers.end();pcl_itr++){
      if((*pcl_itr)->getPdgID()==this->getPdgID())
        return (*pcl_itr)->findProductionMothers();
    }
    return mothers;
}

void TauolaParticle::decay(){

  //Do the decay and set the polarimetric vectors
  TauolaDecay(getSign(),&m_pol_x, &m_pol_y, &m_pol_z, &m_pol_n);
}

void TauolaParticle::addDecayToEventRecord(){

  //Add to decay list used by f_filhep.c
  DecayList::addToEnd(this);
  TauolaWriteDecayToEventRecord(getSign());

  double xmom[4]={0};
  double *pp=xmom;

  if (Tauola::ion[2])
    for(int i=1;1;i++)
    {
      TauolaParticle *x;
      try{ x=DecayList::getParticle(i); }
      catch(std::out_of_range d) {break;}
      if(x->getPdgID()==221){ 

        // Fix 28.04.2011 The pp vector must have boost for eta undone
        TauolaParticle *x_copy = x->clone();
        if(getP(TauolaParticle::Z_AXIS)>0)
          x_copy->boostAlongZ(-getP(),getE());
        else
          x_copy->boostAlongZ(getP(),getE());

        pp[3]=x_copy->getE();
        pp[0]=x_copy->getPx();
        pp[1]=x_copy->getPy();
        pp[2]=x_copy->getPz();
        taueta_(pp,&i);
      } 
    }

  if (Tauola::ion[1])
    for(int i=1;1;i++)
    {
      TauolaParticle *x;
      try{ x=DecayList::getParticle(i); }
      catch(std::out_of_range d) {break;}
      if(x->getPdgID()==310){

        // Fix 28.04.2011 The pp vector must have boost for k0 undone
        TauolaParticle *x_copy = x->clone();
        if(getP(TauolaParticle::Z_AXIS)>0)
          x_copy->boostAlongZ(-getP(),getE());
        else
          x_copy->boostAlongZ(getP(),getE());

        pp[3]=x_copy->getE();
        pp[0]=x_copy->getPx();
        pp[1]=x_copy->getPy();
        pp[2]=x_copy->getPz();
        tauk0s_(pp,&i);
      } 
    }

  if (Tauola::ion[0])
    for(int i=1;1;i++)
    {
      TauolaParticle *x;
      try{ x=DecayList::getParticle(i); }
      catch(std::out_of_range d) {break;}
      if(x->getPdgID()==111){ 

        // Fix 28.04.2011 The pp vector must have boost for pi0 undone
        TauolaParticle *x_copy = x->clone();
        if(getP(TauolaParticle::Z_AXIS)>0)
          x_copy->boostAlongZ(-getP(),getE());
        else
          x_copy->boostAlongZ(getP(),getE());

        pp[3]=x_copy->getE();
        pp[0]=x_copy->getPx();
        pp[1]=x_copy->getPy();
        pp[2]=x_copy->getPz();
        taupi0_(pp,&i);
      } 
    }
  DecayList::clear();

  if(!hasDaughters())
    Log::Fatal("TAUOLA failed. No decay was created",5);
  //  checkMomentumConservation();
  //  decayEndgame(); // vertex shift was wrongly calculated, 
  //                     used 4-momenta should be in lab frame,
  //                     thanks to Sho Iwamoto for debug info.

}


void TauolaParticle::boostDaughtersFromRestFrame(TauolaParticle * tau_momentum){

  if(!hasDaughters()) //if there are no daughters
    return;

  vector<TauolaParticle*> daughters = getDaughters();
  vector<TauolaParticle*>::iterator pcl_itr = daughters.begin();
  
  //get all daughters then rotate and boost them.
  for(;pcl_itr != daughters.end();pcl_itr++){
 
    (*pcl_itr)->boostFromRestFrame(tau_momentum);
    (*pcl_itr)->boostDaughtersFromRestFrame(tau_momentum);
  }
  //checkMomentumConservation();
}

void TauolaParticle::boostDaughtersToRestFrame(TauolaParticle * tau_momentum){

  if(!hasDaughters()) //if there are no daughters
    return;
  // NOTE: Not executed by release examples
  //       because !hasDaughters() is always true
  vector<TauolaParticle*> daughters = getDaughters();
  vector<TauolaParticle*>::iterator pcl_itr = daughters.begin();
  
  //get all daughters then rotate and boost them.
  for(;pcl_itr != daughters.end();pcl_itr++){
 
    (*pcl_itr)->boostToRestFrame(tau_momentum);
    (*pcl_itr)->boostDaughtersToRestFrame(tau_momentum);
  }
  //checkMomentumConservation();
}


void TauolaParticle::boostToRestFrame(TauolaParticle * tau_momentum){

  double theta = tau_momentum->getRotationAngle(Y_AXIS);
  tau_momentum->rotate(Y_AXIS,theta);
  double phi = tau_momentum->getRotationAngle(X_AXIS);
  tau_momentum->rotate(Y_AXIS,-theta);

  //Now rotate coordinates to get boost in Z direction.
  rotate(Y_AXIS,theta);
  rotate(X_AXIS,phi);
  boostAlongZ(-1*tau_momentum->getP(),tau_momentum->getE());
  rotate(X_AXIS,-phi);
  rotate(Y_AXIS,-theta);

}

void TauolaParticle::boostFromRestFrame(TauolaParticle * tau_momentum){
  //get the rotation angles
  //and boost z

  double theta = tau_momentum->getRotationAngle(Y_AXIS);
  tau_momentum->rotate(Y_AXIS,theta);
  double phi = tau_momentum->getRotationAngle(X_AXIS);
  tau_momentum->rotate(Y_AXIS,-theta);

  //Now rotate coordinates to get boost in Z direction.
  rotate(Y_AXIS,theta);
  rotate(X_AXIS,phi);
  boostAlongZ(tau_momentum->getP(),tau_momentum->getE());
  rotate(X_AXIS,-phi);
  rotate(Y_AXIS,-theta);
}

/** Get the angle needed to rotate the 4 momentum vector so that
    the x (y) component disapears. (and the Z component is > 0) */
double TauolaParticle::getRotationAngle(int axis, int second_axis){

  /**if(getP(axis)==0){
    if(getPz()>0)
      return 0; //no rotaion required
    else
      return M_PI;
      }**/
  if(getP(second_axis)==0){
    if(getP(axis)>0)
      return -M_PI/2.0;
    else
      return M_PI/2.0;
  }
  if(getP(second_axis)>0)
    return -atan(getP(axis)/getP(second_axis));
  else
    return M_PI-atan(getP(axis)/getP(second_axis));

}

/** Boost this vector along the Z direction.
    Assume no momentum components in the X or Y directions. */
void TauolaParticle::boostAlongZ(double boost_pz, double boost_e){

  // Boost along the Z axis
  double m=sqrt(boost_e*boost_e-boost_pz*boost_pz);

  double p=getPz();
  double e=getE();

  setPz((boost_e*p + boost_pz*e)/m);
  setE((boost_pz*p + boost_e*e )/m);
}

/** Rotation around an axis X or Y */
void TauolaParticle::rotate(int axis,double theta, int second_axis){
  
  double temp_px=getP(axis);
  double temp_pz=getP(second_axis);
  setP(axis,cos(theta)*temp_px + sin(theta)*temp_pz);
  setP(second_axis,-sin(theta)*temp_px + cos(theta)*temp_pz);
}

void TauolaParticle::rotateDaughters(int axis,double theta, int second_axis){
  if(!hasDaughters()) //if there are no daughters
    return;

  vector<TauolaParticle*> daughters = getDaughters();
  vector<TauolaParticle*>::iterator pcl_itr = daughters.begin();
  
  //get all daughters then rotate and boost them.
  for(;pcl_itr != daughters.end();pcl_itr++){
 
    (*pcl_itr)->rotate(axis,theta,second_axis);
    (*pcl_itr)->rotateDaughters(axis,theta,second_axis);
  }
  //checkMomentumConservation();
}

double TauolaParticle::getMass(){
  double e_sq=getE()*getE();
  double p_sq=getP()*getP();

  if(e_sq>p_sq)
    return sqrt(e_sq-p_sq);
  else
    return -1*sqrt(p_sq-e_sq); //if it's negative
}

double TauolaParticle::getP(){
  return sqrt(getPx()*getPx()+getPy()*getPy()+getPz()*getPz());
}

double TauolaParticle::getP(int axis){
  if(axis==X_AXIS)
    return getPx();

  if(axis==Y_AXIS)
    return getPy();

  if(axis==Z_AXIS)
    return getPz();

  return 0;
}

void TauolaParticle::setP(int axis, double p_component){
  if(axis==X_AXIS)
    setPx(p_component);
  if(axis==Y_AXIS)
    setPy(p_component);
  if(axis==Z_AXIS)
    setPz(p_component);
}

} // namespace Tauolapp
