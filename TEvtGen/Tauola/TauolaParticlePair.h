#ifndef _TauolaParticlePair_h_included_
#define _TauolaParticlePair_h_included_

/**
 * @class TauolaParticlePair
 *
 * @brief Contains two TauolaParticle that are related by 
 * the same mother. Spin correlations are handled here.
 *
 * An object of TauolaParticlePair contains two TauolaParticle
 * that are related by the same mother. Generally this will be
 * a tau+ and tau- or a tau and tau neutrino. For the case of
 * event records that contain multiple instances of the same 
 * particle. eg. tau -> gamma tau or simply tau -> tau. Both
 * the tau from the production vertex, and the final tau before
 * the decay vertex are stored. This allows better handling
 * of spin correlations. The decay is done in the rest frame of 
 * the final tau, where as the spin weight is calculated in the
 * rest frame of the production tau. All spin weights are done
 * in this class. Please refer to the decayTauPairs() method.
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 */


#include <iostream>
#include <vector>
#include <math.h>
#include "TauolaParticle.h"

namespace Tauolapp
{

class TauolaParticlePair{

 public:

  //needed to access m_R matrix and recalculateRij() function.
  friend class Plots;

  /** This constructor takes the TauolaParticle and traverse
      the event structure to find the mother, partner tau or tau 
      neutrino and assosiated final and production versions.
      Once a TauolaParticlePair object has been created in this way
      it is ready to be decayed via decayTauPairs(). */
  TauolaParticlePair(std::vector<TauolaParticle*> &particle_list);

  /** Call the decay method of each 'final' tau. Then calculate
      the spin correlation weight from the particles polarimetric 
      vectors. Decays are accepted or rejected based on the spin 
      weight. Rejected decays are redecayed. */
  void decayTauPair();

  /** Does this pair contain the particle "particle". Note: it only
   checks the "final" particles. */
  bool contains(TauolaParticle * particle);

  /** Return the tau+ particle */
  TauolaParticle * getTauPlus(std::vector<TauolaParticle*> particles);

  /** Return the tau- particle */
  TauolaParticle * getTauMinus(std::vector<TauolaParticle*> particles);

  /** Return the first grandmother of the tau-
      which is an anti-quark or anti-lepton. */
  TauolaParticle * getGrandmotherPlus(std::vector<TauolaParticle*> particles);

  /** Return the first grandmother of the tau-
      which is a quark or lepton. */
  TauolaParticle * getGrandmotherMinus(std::vector<TauolaParticle*> particles);

  /** Print information about the mother and tau pair (at production and final). */
  void print();

  /** Check that the 4 momentum in conserved at the verticle of
      each decayed tau. */
  void checkMomentumConservation();

 private:

  /** Default constructor is private, so that only friend class can use it. */
  TauolaParticlePair() {}

  /** Store born variables in Tauola class, so the user can retrieve
      them using Tauola::getBornKinematics. */
  static void setBornKinematics(int  incoming_pdg_id, int  outgoing_pdg_id, double  invariant_mass_squared, double  cosTheta);

  /** Pointers to taus (or tau and neutrino) as they
      are before being decayed. */
  std::vector<TauolaParticle*> m_final_particles;
  
  /** Pointers to taus (or tau and neutrino) as they
      are after production. */
  std::vector<TauolaParticle*> m_production_particles;
  
  /** Pointer to mothers of the tau pair. */  
  TauolaParticle* m_mother;

  /** Is there an entry in the event record for the tau pair's mother? */  
  bool m_mother_exists;
  
  /** vector of pointers to the taus grandparents */
  std::vector<TauolaParticle*> m_grandmothers;

  /** If SANC tables are present, use them to recalculate the matrix Rij. */
  void recalculateRij(int incoming_pdg_id, int outgoing_pdg_id, double invariant_mass_squared, double cosTheta);

  /** Rotate the whole system using the given angle theta. */
  void rotateSystem(vector<TauolaParticle *> grandmothers,
                    vector<TauolaParticle *> taus,
                    double theta,
                    int axis, 
                    int axis2=TauolaParticle::Z_AXIS);


  /** Boost the outgoing tau and partner and the incoming grandparents of
      the tau to the mothers rest frame. The mother is not boosted.
      The axis are rotated so that the particle given by "z_axis_particle" is aligned
      on the z-axis. If "alignment" is -1 is will be aligned in the negative z direction.
      otherwise it is aligned in the positive direction. rotaion_angle(1-3) are
      returned to allow reversal of the transformation (through the method 
      boostFromMotherToLabFrame).*/
  void boostFromLabToTauPairFrame(double * rotation_angle1, 
                                  double * rotation_angle2,
                                  double * rotation_angle3,
                                  TauolaParticle * mother,
                                  vector<TauolaParticle *> grandmothers,
                                  vector<TauolaParticle *> taus);

  /** Reverses the transformation of boostFromLabToMothersFrame. **/
  void boostFromTauPairToLabFrame(double rotation_angle1, 
                                  double rotation_angle2,
                                  double rotation_angle3,
                                  TauolaParticle * mother,
                                  vector<TauolaParticle *> grandmothers,
                                  vector<TauolaParticle *> taus);
    
  /** The density matric m_R is filled based on the mothers type and kinematics
      of the event in the mothers rest frame. */
  void initializeDensityMatrix();

  /** create a particle which m_mother points to. This is based on the
      daughters 4-momentum and particle type. A Z or W is assumed if the
      configuration of taus and neutrinos is correct. This particle is not
      written into the event record, but it used by the fillDenistyMatrix
      method for spin correlations */
  TauolaParticle * makeTemporaryMother(vector<TauolaParticle *> taus);

  /**Needs to be changed*/
  //void ANGULU(int *IDE, int *IDF, double *SVAR, double *COSTHE);
  /**Needs to be changed*/
  double getZPolarization(int *incoming_pdg_id,
                          int *outgoing_pdg_id, 
                          double *invMass,
                          double *cosTheta);

  /** Private function, calculates virtuality between two particles. */
  double getVirtuality(TauolaParticle * p1, TauolaParticle*p2, bool flip);

  /** Add particle to beam. */
  void addToBeam(TauolaParticle * pcle,
                 std::vector<TauolaParticle*> * candidates_same,
                 std::vector<TauolaParticle*> * candidates_opp);
    

 /** frames in which it is defined are fixed by the methods 
      boostFromLabToMotherFrame and boostFromMotherToLabFrame.
      Modification to m_R and boostFrom/ToMotherFrame must be done
      coherently. */
  double m_R[4][4]; //density matrix
};

//Temporary
//Pz is still calculated using the FORTRAN routine in tauola_extra.f
//This should be migrated to C++ at some stage.
extern "C" {
  extern double plzap0_(int *incoming_pdg_id,
                        int *outgoing_pdg_id, 
                        double *invMass,
                        double *cosTheta);
}

} // namespace Tauolapp
#endif  

