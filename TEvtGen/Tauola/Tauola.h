#ifndef _Tauola_h_included_
#define _Tauola_h_included_

/** 
 * @class Tauola
 *
 * @brief Controls the configuration, initialization of Tauola.
 *
 * The Tauola class provides a wrapper to the TAUOLA common block
 * responsible for setting up TAUOLA. They should then configure Tauola
 * via the set method and then call initialize(). This is a static class.
 *
 * @author Nadia Davidson
 * @date 16th June 2008
 */

#include <iostream>
#include "TauolaParticle.h"
#include "f_InitTauola.h"
#include "f_Variables.h"

namespace Tauolapp
{

class TauolaEvent;
class TauolaParticle;

class Tauola{

 public:

  static const int NS1=100,NS2=100,NS3=100,NCOS=21;

  /** Units */
  static enum MomentumUnits { DEFAULT_MOMENTUM=-1, MEV, GEV } momentumUnit;
  static enum LengthUnits   { DEFAULT_LENGTH  =-1, MM , CM  } lengthUnit;

  /** Set output units (default Tauola::GEV and Tauola::MM). */
  static void   setUnits(MomentumUnits m,LengthUnits l);

  /** Set tau lifetime (in mm). */
  static void   setTauLifetime(double t);

  /** Decay Modes */
  enum { All=0, ElectronMode, MuonMode, PionMode,
         RhoMode, A1Mode, KMode, KStarMode };

  /** Structure for switching the computation of spin correlation.
      By default all spin correlations are turned on. */
  struct Particles
  {
     bool GAMMA,
          Z0,
          HIGGS,
          HIGGS_H,
          HIGGS_A,
          HIGGS_PLUS,
          HIGGS_MINUS,
          W_PLUS,
          W_MINUS;
     void setAll(bool flag) { GAMMA=Z0=HIGGS=HIGGS_H=HIGGS_A=HIGGS_PLUS=HIGGS_MINUS=W_PLUS=W_MINUS=flag; }
  } static spin_correlation;

   /** Initalize Tauola with the parameters previously set via the
       setter methods */
   static void initialize();

   /** DEPRECATED: Use 'initialize' instead. */
   static void initialise();

   /** Change currents used by Tauola.
       mode = 0 (default) - use CLEO currents
       mode = 1 use RChL currents for 3pi and Belle currents for 2pi */
   static void setNewCurrents(int mode);

   /** Substitute build-in generator with external one */
   static void setRandomGenerator( double (*gen)() );

   static void setRedefineTauMinus( void (*fun)(TauolaParticle *) );
   static void setRedefineTauPlus ( void (*fun)(TauolaParticle *) );

   /** Tau gun. Takes one particle that's already inside an event record and produces it's decay.
       The tau provided may be undecayed first, or left intact if it already has daughters.
       If the polarization three-vector is provided it will be used to construct m_R matrix. */
   static void decayOne(TauolaParticle *tau, bool undecay=false, double polx=0,double poly=0, double polz=0);

   /** Checks if we are using decayOne() */
   static bool isUsingDecayOne();

   /** Checks if we are using boost routine for decayOne */
   static bool isUsingDecayOneBoost();

   /** Set boost routine for decayOne(). Refer to documentation for more details. */
   static void setBoostRoutine( void (*boost)(TauolaParticle*, TauolaParticle *) );

   /** Execute boost routine for decayOne() */
   static void decayOneBoost(TauolaParticle *mother, TauolaParticle *target);

   /** Return polarization vector used by decayOne() */
   static const double* getDecayOnePolarization();

   /** Set the pdg id of the particle to decay (should be 15 or -15) */
   static void setDecayingParticle(int pdg_id);

   /** Return the pdg id of the particle to decay */
   static int getDecayingParticle();
   
   /** Set the decay mode of all particle with pdg id the same
       as the one given in setDecayingParticle(). firstDecayMode=0 
       is default and allows all decay modes. */
   static void setSameParticleDecayMode(int firstDecayMode);
   
   /** Set the decay mode of all particle with opposite charge
       to the one given in setDecayingParticle(). secondDecayMode=0 
       is default and allows all decay modes. */
   static void setOppositeParticleDecayMode(int secondDecayMode);

   /** Switch for bremssthahlung in leptonic tau decays */
   static void setRadiation(bool rad);

   /** Cut-Off parameter of radition. Above that value photon is explicitly generated */
   static void setRadiationCutOff(double rad_cut_off);

   /** Initialization of some constants related to QED corrections.
       Variable iniphy_param is at present dummy. It is prepared to be transmitted
       to some old style production code and is kept for backward compatibility */
   static void setInitializePhy(double iniphy);

   /** DEPRECATED: Use 'setInitializePhy' instead. */
   static void setInitialisePhy(double iniphy);

   /** Set branching fraction for i-th channel. Can be reused several times during the run. */
   static void setTauBr(int i, double value);

   static void setTaukle(double bra1, double brk0, double brk0b, double brks);

   static double getHiggsScalarPseudoscalarMixingAngle();

   /** set the mixing angle. coupling: tau~(cos(phi)+isin(phi)gamma5)tau */
   static void setHiggsScalarPseudoscalarMixingAngle(double angle); 

   /** Get mass of the tau used by interface. */
   static double getTauMass();

   /** Modify Higgs Scalar-Pseudoscalar PDG id (default is 35). */
   static void setHiggsScalarPseudoscalarPDG(int pdg_id);

   /** Get Higgs Scalar-Pseudoscalar PDG id. */
   static int getHiggsScalarPseudoscalarPDG();

   static int getHelPlus();

   static int getHelMinus();

   static double getEWwt();

   static double getEWwt0();

   static void setEWwt(double wt, double wt0);

   static void setHelicities(int Minus, int Plus);

   static void setEtaK0sPi(int eta, int k, int pi);

   static void getBornKinematics(int *incoming_pdg_id, int *outgoing_pdg_id, double *invariant_mass_squared,double *cosTheta);

   static void summary();

public:

   static double table11A[NS1][NCOS][4][4],table1A[NS1][NCOS][4][4],table2A[NS1][NCOS][4][4];
   static double wtable11A[NS1][NCOS],wtable1A[NS1][NCOS],wtable2A[NS1][NCOS];
   static double w0table11A[NS1][NCOS],w0table1A[NS1][NCOS],w0table2A[NS1][NCOS];

   static double table11B[NS2][NCOS][4][4],table1B[NS2][NCOS][4][4],table2B[NS2][NCOS][4][4];
   static double wtable11B[NS2][NCOS],wtable1B[NS2][NCOS],wtable2B[NS2][NCOS];
   static double w0table11B[NS2][NCOS],w0table1B[NS2][NCOS],w0table2B[NS2][NCOS];

   static double table11C[NS3][NCOS][4][4],table1C[NS3][NCOS][4][4],table2C[NS3][NCOS][4][4];
   static double wtable11C[NS3][NCOS],wtable1C[NS3][NCOS],wtable2C[NS3][NCOS];
   static double w0table11C[NS3][NCOS],w0table1C[NS3][NCOS],w0table2C[NS3][NCOS];
   static double sminA,smaxA,sminB,smaxB,sminC,smaxC;

   static int ion[3];

   // c*tau in milimeters, survival probablility  P(t)=exp(-t/lifetime) 
   static double tau_lifetime;
   static double momentum_conservation_threshold;

   //born kinematic variables
   static int buf_incoming_pdg_id,  buf_outgoing_pdg_id;
   static double  buf_invariant_mass_squared,  buf_cosTheta;
   static double  buf_R[4][4]; //density matrix

   //pointer to random generator function
   static double (*randomDouble)();

   static void (*redefineTauPlusProperties)(TauolaParticle *);
   static void (*redefineTauMinusProperties)(TauolaParticle *);

 private:

  /** Calculate the charge  of particle  with code 'idhep'. 
      The code  of the  particle  is  defined by the Particle Data
      Group in Phys. Lett. B204 (1988) 1.
      NOTE: Code taken from Photos++, file: PhotosUtilities.cxx, function: PHOCHA */
  static double particleCharge(int idhep);
  
  /** Fill 'array' indices from 'beg' to 'end' with 'value' */
  static void fill_val(int beg, int end, double* array, double value);
  
  /** Default generator used in Tauola */
  static double defaultRandomGenerator();
  static void   defaultRedPlus(TauolaParticle *);
  static void   defaultRedMinus(TauolaParticle *);

  /** Are we using decayOne() ? */
  static bool   m_is_using_decay_one;
  /** decayOne() polarization vector */
  static double m_decay_one_polarization[3];
  /** Boost routine used by decayOne() */
  static void (*m_decay_one_boost_routine)(TauolaParticle*,TauolaParticle*);

  static int m_pdg_id;
  static int m_firstDecayMode; 
  static int m_secondDecayMode;
  static bool m_rad;
  static double m_rad_cut_off;
  static double m_iniphy;
  static double m_higgs_scalar_pseudoscalar_mix;
  static int m_higgs_scalar_pseudoscalar_pdg;
  static double m_wtEW;
  static double m_wtEW0;
  static int m_helPlus;
  static int m_helMinus;
};

} // namespace Tauolapp
#endif  

