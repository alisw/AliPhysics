#ifndef _TauolaParticle_h_included_
#define _TauolaParticle_h_included_

/**
 * @class TauolaParticle
 *
 * @brief Abstract base class for particle in the event. This class also
 * handles boosting.
 *
 * TauolaParticle is a Tauola representation of a particle. It has virtual
 * getter and setter methods that need to be implemented by a derived class. 
 * An example of this is TauolaHepMCParticle. In this way it provides an
 * interface to the information in the Event Record.
 *
 * The class is also responsible for decays and contains the polarimetric 
 * vector returned from tauola. All boosting is also done here.
 *
 * @author Nadia Davidson
 * @date 16 June 2008
 */

#include <iostream>
#include <math.h>
#include <vector>

#include "DecayList.h"
#include "Tauola.h"
#include "f_Decay.h"

namespace Tauolapp
{

class TauolaParticle{

 public:

  virtual ~TauolaParticle(){};

  /** The same sign as decaying particle pdg ID code 
      given to Tauola object (only meaningful for taus). */
  static const int SAME_SIGN=1;

  /** The opposite sign to decaying particle pdg ID code 
      given to Tauola object (only meaningful for taus). */
  static const int OPPOSITE_SIGN=2;

  /** Sign type is not applicable for this particle
      (probably it's not a tau). */
  static const int NA_SIGN=3;

  /** Stable particle status */
  static const int STABLE=1;

  /** Decayed particle status */
  static const int DECAYED=2;

  /** History particle status */
  static const int HISTORY=3;  

  /** X Axis */
  static const int X_AXIS=1;

  /** Y Axis */
  static const int Y_AXIS=2;

  /** Z Axis */
  static const int Z_AXIS=3;

  /** Z0 particle */
  static const int Z0 = 23;

  /** H particle */
  static const int HIGGS = 25;

  /** A0 particle */
  static const int HIGGS_A = 36;

  /** H+ particle */
  static const int HIGGS_PLUS = 37;

  /** H- particle */
  static const int HIGGS_MINUS = -37;


  /** W+ particle */
  static const int W_PLUS = 24;

  /** W- particle */
  static const int W_MINUS = -24;

  /** photon */
  static const int GAMMA = 22;

  /** tau+ particle */
  static const int TAU_PLUS = -15;

  /** tau- particle */
  static const int TAU_MINUS = 15;

  /** tau neutrino particle */
  static const int TAU_NEUTRINO = 16;

  /** tau antineutrino particle */
  static const int TAU_ANTINEUTRINO = -16;


  /** muon+ particle */
  static const int MUON_PLUS = -13;

  /** muon- particle */
  static const int MUON_MINUS = 13;

  /** muon neutrino particle */
  static const int MUON_NEUTRINO = 14;

  /** muon antineutrino particle */
  static const int MUON_ANTINEUTRINO = -14;


  /** e+ particle */
  static const int POSITRON = -11;

  /** e- particle */
  static const int ELECTRON = 11;

  /** e neutrino particle */
  static const int ELECTRON_NEUTRINO = 12;

  /** e antineutrino particle */
  static const int ELECTRON_ANTINEUTRINO = -12;

  /** up quark */
  static const int UP = 2;

  /** anti-up quark */
  static const int ANTIUP = -2;

  /** down quark */
  static const int DOWN = 1;

  /** anti-down quark */
  static const int ANTIDOWN = -1;

  static const int GLUON = 21;
  static const int CHARM = 4;
  static const int TOP = 6;
  static const int STRANGE = 3;
  static const int BOTTOM = 5;

  /** All other particle types*/
  static const int OTHER = 0;

  /** Create a new particle with the same properties as this one.
      Mothers and daughters will not be linked. */
  TauolaParticle * clone();

  /** Get the angle between this particle and another particle */
  double getAngle(TauolaParticle *);

  /** Add the 4 momentum of another particle to this particle */
  void add(TauolaParticle *);

  /** Subtract the 4 momentum of another particle from this particle */
  void subtract(TauolaParticle *);

  /** Decay the particle. This calls the decay methods in the
      interface to (FORTRAN) tauola. */
  void decay();

  /** Invokes TAUOLA FORTRAN routine DEKAY retrieving the daughters of
      decayed tau. */
  void addDecayToEventRecord();

  /** Get whether this particle has the same or opposite sign
      to the pdg code given to the Tauola object. (Only relevant
      for tau particles).*/  
  int  getSign();

  /** Get the polarimetric vector of this particle in the direction X. 
      (Only relevant for tau particles).*/
  double getPolarimetricX();

  /** Get the polarimetric vector of this particle in the direction Y. 
      (Only relevant for tau particles).*/
  double getPolarimetricY();

  /** Get the polarimetric vector of this particle in the direction Z. 
      (Only relevant for tau particles).*/
  double getPolarimetricZ(); 

  /** Return whether the particle has any chidren */
  bool hasDaughters();

  /** Traverse the event structure and find the final version
      of this particle which does not have a particle of it's own type
      as it's daughter. eg. Generally the final stable copy */ 
  TauolaParticle * findLastSelf();

  /** Traverse the event structure and find the first set of mothers 
      which are not of the same type as this particle. */
  std::vector<TauolaParticle *> findProductionMothers();

  /** Transform this particles four momentum from the lab frome
      into the rest frame of the paramter TauolaParticle. **/   
  void boostToRestFrame(TauolaParticle * boost);

  /** Transform the four momentum of all the daughters recursively 
      into the frame of the "particle" TauolaParticle. **/   
  void boostDaughtersToRestFrame(TauolaParticle * boost);


  /** Transform this particles four momentum from the rest frame of
      the paramter TauolaParticle, back into the lab frame. **/   
  void boostFromRestFrame(TauolaParticle * boost);

  void boostDaughtersFromRestFrame(TauolaParticle * boost);

  /** Do a Lorenz transformation along the Z axis. */ 
  void boostAlongZ(double pz, double e);

  /** rotate this particles 4-momentum by an angle phi from
   the axisis "axis" towards the axis "second_axis". */
  void rotate(int axis, double phi, int second_axis=Z_AXIS);

  void rotateDaughters(int axis, double phi, int second_axis=Z_AXIS);

  /** Returns the angle around the axis "axis" needed to rotate 
      the four momenum is such a way that the non-Z component 
      disappears and Z>0. This is used to rotating the coordinate 
      system into a frame with only a Z component before calling 
      boostAlongZ().*/
  double getRotationAngle(int axis, int second_axis=Z_AXIS);

  /** Get scalar momentum */
  double getP();

  /** Get momentum component in the direction of "axis" (x,y,z) */
  double getP(int axis);

  /** Set momentum component in the direction of "axis" (x,y,z) */
  void  setP(int axis, double p_component);

  /** Get the invariant mass from the four momentum*/
  double getMass();



  /********************************************** 
      Beginning of virtual methods 

  ********************************************/

  /** remove the ougoing branch from this particles and reset its status to stable */
  virtual void undecay(){};

  /** check that the 4 momentum in conserved at the vertices producing
      and ending this particle */
  virtual void checkMomentumConservation(){};

  /** Optional. Modify particle or decay tree if needed. */
  virtual void decayEndgame(){};

  /** Returns the px component of the four vector*/
  virtual double getPx()=0;

  /** Returns the py component of the four vector */
  virtual double getPy()=0;

  /** Returns the pz component of the four vector */
  virtual double getPz()=0;

  /** Returns the energy component of the four vector */
  virtual double getE()=0;

  /** Set the px component of the four vector */
  virtual void setPx( double px )=0;

  /** Set the px component of the four vector */
  virtual void setPy( double py )=0;

  /** Set the pz component of the four vector */
  virtual void setPz( double pz )=0;

  /** Set the energy component of the four vector */
  virtual void setE( double e )=0;

  /** Set the mothers of this particle via a vector of TauolaParticle */
  virtual void setMothers(std::vector<TauolaParticle*> mothers)=0;

  /** Set the daughters of this particle via a vector of TauolaParticle */  
  virtual void setDaughters(std::vector<TauolaParticle*> daughters)=0;

  /** Returns the mothers of this particle via a vector of TauolaParticle */
  virtual std::vector<TauolaParticle*> getMothers()=0;

  /** Returns the daughters of this particle via a vector of TauolaParticle */
  virtual std::vector<TauolaParticle*> getDaughters()=0;

  /** Set the PDG ID code of this particle */
  virtual void setPdgID(int pdg_id)=0;

  /** Set the mass of this particle */
  virtual void setMass(double mass)=0;

  /** Set the status of this particle */
  virtual void setStatus(int status)=0;

  /** Get the PDG ID code of this particle */
  virtual int getPdgID()=0;

  /** Get the status of this particle */
  virtual int getStatus()=0;

  /** Get the barcode of this particle */
  virtual int getBarcode()=0;

  /** Create a new particle of the same type, with the given
      properties. The new particle bares no relations to this
      particle, but it provides a way of creating a intance of
      the derived class. eg. createNewParticle() is used inside
      filhep_() so that an eg. TauolaHepMCParticle is created without
      the method having explicit knowledge of the TauolaHepMCParticle 
      class */
  virtual TauolaParticle * createNewParticle(int pdg_id, int status, 
                                             double mass, double px,
                                             double py, double pz,
                                             double e)=0;

  /** Print some information about this particle to standard output */
  virtual void print()=0;

 private:

  /** The polarimetric vector of this particle in the direction X. 
      (Only relevant for tau particles). */
  double m_pol_x;

  /** The polarimetric vector of this particle in the direction Y. 
      (Only relevant for tau particles). */
  double m_pol_y;

  /** The polarimetric vector of this particle in the direction Z. 
      (Only relevant for tau particles). */
  double m_pol_z;

  /** Fourth component of the polarimetric vector. Should be the
      normalisation (1). (Only relevant for tau particles). */
  double m_pol_n;
};

} // namespace Tauolapp
#endif  

