// HIUserHooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the definition of the HIUserHooks class and a
// set of other classes used inside Pythia8 to model collisions
// involving heavy ions.
// Nucleon: represents a proton or a neutron inside a necleus.
// SubCollision: a collision between a projectile and a target Nucleon.
// NucleusModel: models the Nucleon distribution in a nucleus.
// WoodsSaxonModel: NucleusModel implementing a simple Woods-Saxon.
// GLISSANDOModel: NucleusModel implementing the GLISSANDO prescription.
// ImpactParameterGenerator: distributes nuclei in impact parameter space.
// SubCollisionModel: Models the collision probabilities of nucleons.
// NaiveSubCollisionModel: A very simple SubCollisionModel.
// DoubleStrikman: A more advanced SubCollisionModel.
// EventInfo: stores full nucleon-nucleon events with corresponding Info.
// HIInfo: info about a Heavy Ion run and its produced events.
// HIUserHooks: User hooks for HeavyIons models.

#ifndef Pythia8_HIUserHooks_H
#define Pythia8_HIUserHooks_H

#include "Pythia8/Pythia.h"
#include <limits>

namespace Pythia8 {

class Pythia;

/// Forward declaration.
class EventInfo;

//==========================================================================

// The HIUnits namespace defines the unitsystem used by the heavy ion
// machinery in Pythia8. In particular all lengths are in femtometer
// and cross sections are therefore in squared femtometer.
namespace HIUnits {

// Lengths
const double femtometer = 1.0;
const double millimeter = 1.0e12;

// Cross sections
const double femtometer2 = 1.0;
const double millibarn = 0.1;
const double nanobarn = 1.0e-7;

}

using namespace HIUnits;

//==========================================================================

/// The Nucleon class represent a nucleon in a nucleus. It has an id
/// number (proton or neutron) an impact parameter position (absolute
/// and relative to the nucleus center), a status and a state to be
/// defined and used by a SubCollisionModel.

class Nucleon {

  friend class SubCollisionModel;

public:

  /// Enum for specifying the status of a nucleon.
  enum Status {
    UNWOUNDED = 0,  ///< The nucleon is not wounded.
    ELASTIC = 1,    ///< The nucleon is elastically scattered.
    DIFF = 2,       ///< The nucleon is diffractively wounded.
    ABS = 3         ///< The nucleon is absorptively wounded.
  };

  /// The state of a nucleon is a general vector of doubles.
  typedef vector<double> State;

  /// The constuctor takes a particle id and a position in impact
  /// parameter relative to the nucleus center as arguments.
  Nucleon(int idIn = 0, int indexIn = 0, const Vec4 & pos = Vec4())
    : idSave(idIn), indexSave(indexIn), nPosSave(pos), bPosSave(pos),
      statusSave(UNWOUNDED),eventp(0), isDone(0) {}

  /// Accessor functions:

  /// The nucleon type.
  int id() const { return idSave; }

  /// The nucleon type.
  int index() const { return indexSave; }

  /// The position of this nucleon relative to the nucleus center.
  const Vec4 & nPos() const { return nPosSave; }

  /// The absolute position in impact parameter space.
  const Vec4 & bPos() const { return bPosSave; }

  /// Shift the absolute position in impact parameter space.
  const Vec4 & bShift(const Vec4 & bvec) { return bPosSave += bvec; }

  /// The status.
  Status status() const { return statusSave; }

  /// Check if nucleon has been assigned.
  bool done() const { return isDone; }

  /// The event this nucleon is assigned to.
  EventInfo * event() const { return eventp; }

  /// The physical state of the incoming nucleon.
  const State & state() const { return stateSave; }

  /// Return an alternative state.
  const State & altState(int i = 0) {
    static State nullstate;
    return i < int(altStatesSave.size())? altStatesSave[i]: nullstate;
  }

  /// Manipulating functions:

  /// Set the status.
  void status(Status s) { statusSave = s; }

  /// Set the physical state.
  void state(State s) { stateSave = s; }

  /// Add an alternative state.
  void addAltState(State s) { altStatesSave.push_back(s); }

  /// Select an event for this nucleon.
  void select(EventInfo & evp, Status s) {
    eventp = &evp;
    isDone = true;
    status(s);
  }

  /// Select this nucleon to be assigned to an event.
  void select() { isDone = true; }

  /// Print out debugging information.
  void debug();

private:

  /// The type of nucleon.
  int idSave;

  /// The index of this nucleon.
  int indexSave;

  /// The position in impact parameter relative to the nucleus center.
  Vec4 nPosSave;

  /// The absolute position in impact parameter.
  Vec4 bPosSave;

  /// The status.
  Status statusSave;

  /// The state of this nucleon.
  State stateSave;

  /// Alternative states to be used to understand fluctuations in the
  /// state of this nucleon.
  vector<State> altStatesSave;

  /// Pointer to the even this nucleon ends up in.
  EventInfo * eventp;

  /// True if this nuclein has been assigned to an event.
  bool isDone;

  /// Reset the states and status.
  void reset() {
    statusSave = UNWOUNDED;
    altStatesSave.clear();
    bPosSave = nPosSave;
    isDone = false;
    eventp = 0;
  }

};

//==========================================================================

/// SubCollision represents a possible collision between a projectile
/// and a target Nucleon.

class SubCollision {

public:

  /// This defines the type of a bunary nucleon collison.
  enum Type {
    NONE,       ///< This is not a collision.
    ELASTIC,    ///< This is an elastic scattering
    SDEP,       ///< The projectile is diffractively excited.
    SDET,       ///< The target is diffractively excited.
    DDE,        ///< Both projectile and target are diffractively excited.
    CDE,        ///< Both excited but with central diffraction.
    ABS         ///< This is an absorptive (non-diffractive) collision.
  };

  SubCollision(Nucleon & projIn, Nucleon & targIn,
               double bIn, double bpIn, Type typeIn)
    : proj(&projIn), targ(&targIn), b(bIn), bp(bpIn), type(typeIn) {}

  SubCollision()
    : proj(0), targ(0), b(0.0), bp(0.0), type(NONE) {}

  // Used to order sub-collisions in a set.
  bool operator< (const SubCollision & s) const { return b < s.b; }

  // Return 0 if neither proj or target are neutrons, 1 if target is
  // neutron, 2 if projectile is neutron, and 3 if both are neutrons.
  int nucleons() const {
    return ( abs(targ->id()) == 2112? 1: 0 ) +
           ( abs(proj->id()) == 2112? 2: 0 );
  }

  // The projectile nucleon.
  Nucleon * proj;

  /// The target nucleon.
  Nucleon * targ;

  /// The impact parameter distance between the nucleons in femtometer.
  double b;

  /// The impact parameter distance between the nucleons scaled like
  /// in Pythia to have unit average for non-diffractive collisions.
  double bp;

  /// The type of collison.
  mutable Type type;

};

//==========================================================================

/// This class generates the impact parameter distribution of nucleons
/// in a nucleus.

class NucleusModel {

public:

  /// Default constructor giving the nucleis id and an optional
  /// raduis (in femtometer).
  NucleusModel()
    : idSave(0), ISave(0), ASave(0), ZSave(0), LSave(0), RSave(0.0),
      settingsPtr(0), particleDataPtr(0), rndPtr(0) {}

  /// Virtual destructor.
  virtual ~NucleusModel() {}

  /// Init method.
  void initPtr(int idIn, Settings & settingsIn,
               ParticleData & particleDataIn, Rndm & rndIn);
  virtual bool init();

  virtual Particle produceIon(bool istarg);

  /// Generate a vector of nucleons according to the implemented model
  /// for a nucleus given by the PDG number.
  virtual vector<Nucleon> generate() const = 0;

  /// Accessor functions.
  int id() const { return idSave; }
  int I() const { return ISave; }
  int A() const { return ASave; }
  int Z() const { return ZSave; }
  int L() const { return LSave; }
  double R() const { return RSave; }

protected:

  /// The nucleus.
  int idSave;

  /// Cache information about the nucleus.
  int ISave, ASave, ZSave, LSave;

  /// The estimate of the nucleus radius.
  double RSave;

  /// Pointers to useful objects.
  Settings * settingsPtr;
  ParticleData * particleDataPtr;
  Rndm * rndPtr;

};

//==========================================================================

/// A general Woods-Saxon distributed nucleus.

class WoodsSaxonModel: public NucleusModel {

public:

  /// The default constructor needs a nucleus id, a radius, R, and a
  /// "skin width", a (both length in femtometers).
  WoodsSaxonModel(): aSave(0.0), intlo(0.0),
                    inthi0(0.0), inthi1(0.0), inthi2(0.0) {}

  /// Accessor functions:
  double a() const { return aSave; }

protected:

  /// Generate the position of a single nucleon. (The time component
  /// is always zero).
  Vec4 generateNucleon() const;

  /// Virtual destructor.
  virtual ~WoodsSaxonModel() {}

  /// Setup the generation with a given nucleus radius, R, and a "skin
  /// width", a (both length in femtometers).
  virtual bool init() {
    intlo = R()*R()*R()/3.0;
    inthi0 = a()*R()*R();
    inthi1 = 2.0*a()*a()*R();
    inthi2 = 2.0*a()*a()*a();
    return NucleusModel::init();
  }

protected:

  /// The nucleus radius, skin depth parameter, and hard core nucleon
  /// radius..
  double aSave;

private:

  /// Cashed integrals over the different parts of the over estimating
  /// functions.
  double intlo, inthi0, inthi1, inthi2;

};


//==========================================================================

/// The GLISSANDOModel has a specific parameteraization of a
/// Wood-Saxon potential for A>16 and is described in asXiv:1310.5475
/// [nucl-th].

class GLISSANDOModel: public WoodsSaxonModel {

public:

  /// Default constructor.
  GLISSANDOModel(): RhSave(0.0), gaussHardCore(false) {}

  /// Virtual destructor.
  virtual ~GLISSANDOModel() {}

  /// Initialize.
  bool init();

  /// Generate a vector of nucleons according to the implemented model
  /// for a nucleus given by the PDG number.
  virtual vector<Nucleon> generate() const;

  /// Accessor functions.
  double Rh() const { return RhSave; }

  double RhGauss() const { return RhSave*abs(rndPtr->gauss()); };

private:

  /// The hard core radius;
  double RhSave;

  /// Option to use a Gaussian hard core instead of a sharp one.
  bool gaussHardCore;

};

//==========================================================================

/// Forward Declarations
class SubCollisionModel;
class NucleusModel;

/// ImpactParameterGenerator is able to generate a specific impact
/// parameter together with a weight such that aweighted average over
/// any quantity X(b) corresponds to the infinite integral over d^2b
/// X(b). This base class gives a Gaussian profile, d^2b exp(-b^2/2w^2).

class ImpactParameterGenerator {

public:

  /// The default constructor takes a gneral width (in femtometers) as
  /// argument.
  ImpactParameterGenerator()
    : widthSave(0.0), collPtr(0), projPtr(0), targPtr(0),
      settingsPtr(0), rndPtr(0) {}

   /// Virtual destructor.
  virtual ~ImpactParameterGenerator() {}

  /// Virtual init method.
  virtual bool init();
  void initPtr(SubCollisionModel & collIn,
               NucleusModel & projIn,
               NucleusModel & targIn,
               Settings & settingsIn,
               Rndm & rndIn);

  /// Return a new impact parameter and set the corresponding weight
  /// provided.
  virtual Vec4 generate(double & weight) const;

  /// Set the width (in femtometers).
  void width(double widthIn) { widthSave = widthIn; }

  /// Get the width.
  double width() const { return widthSave; }

private:

  /// The width of a distribution.
  double widthSave;

protected:

  /// Info from the controlling HeavyIons object
  SubCollisionModel * collPtr;
  NucleusModel * projPtr;
  NucleusModel * targPtr;
  Settings * settingsPtr;
  Rndm * rndPtr;

};


//==========================================================================

/// The SubCollisionModel is is able to model the collision between two
/// nucleons to tell which type of collision has occurred. The model
/// may manipulate the corresponing state of the nucleons.

class SubCollisionModel {

public:

  /// Internal class to report cross section estimates.
  struct SigEst {
    /// The cross sections (tot, nd, dd, sdp, sdt, cd, el, bslope).
    vector<double> sig;

    /// The extimated error (squared)
    vector<double> dsig2;

    /// Which cross sections were actually fitted
    vector<bool> fsig;

    /// The estimate of the average (and squared error) impact
    /// parameter for inelastic non-diffractive collisions.
    double avNDb, davNDb2;

    /// Constructor for zeros.
    SigEst(): sig(8, 0.0), dsig2(8, 0.0), fsig(8, false),
              avNDb(0.0), davNDb2(0.0) {}

  };

public:

  /// The default constructor is empty.
  SubCollisionModel(): sigTarg(8, 0.0), sigErr(8, 0.05), NInt(100000),
    NGen(20), NPop(20), sigFuzz(0.2), fitPrint(true), avNDb(1.0*femtometer),
    projPtr(), targPtr(), sigTotPtr(), settingsPtr(), infoPtr(), rndPtr() {}

  /// Virtual destructor,
  virtual ~SubCollisionModel() {}

  /// Virtual init method.
  virtual bool init();

  void initPtr(NucleusModel & projIn, NucleusModel & targIn,
               SigmaTotal & sigTotIn, Settings & settingsIn,
               Info & infoIn, Rndm & rndIn) {
    projPtr = &projIn;
    targPtr = &targIn;
    sigTotPtr = &sigTotIn;
    settingsPtr = &settingsIn;
    infoPtr = &infoIn;
    rndPtr = & rndIn;
  }

  /// Take two vectors of Nucleons and an impact parameter vector and
  /// produce the corrsponding sub-collisions. Note that states of the
  /// nucleons may be changed. The function in this abstract base
  /// class will reset the nucleon states for convenience. The
  /// sub-collisions are ordered in the impact parameter distance
  /// between the nucleons. The T-variable will be set to the summed
  /// elastic amplityde.
  virtual multiset<SubCollision> getCollisions(vector<Nucleon> & proj,
                                               vector<Nucleon> & targ,
                                               const Vec4 & bvec,
                                               double & T) = 0;

  /// Access the nucleon-nucleon cross sections assumed
  /// for this model.

  /// The total cross section.
  double sigTot() const {
    return sigTarg[0];
  }

  /// The total cross section.
  double sigEl() const { return sigTarg[6]; }

  /// The central diffractive excitation cross section.
  double sigCDE() const { return sigTarg[5]; }

  /// The single diffractive excitation cross section (both sides summed).
  double sigSDE() const { return sigTarg[3] + sigTarg[4]; }

  /// The single diffractive excitation cross section (excited projectile).
  double sigSDEP() const { return sigTarg[3]; }

  /// The single diffractive excitation cross section (excited target).
  double sigSDET() const { return sigTarg[4]; }

  /// The double diffractive excitation cross section.
  double sigDDE() const { return sigTarg[2]; }

  /// The non-diffractive (absorptive) cross section.
  double sigND() const { return sigTarg[1]; }

  /// The elastic b-slope parameter.
  double bSlope() const { return sigTarg[7]; }

  /// Calculate the cross sections for the given set of parameters.
  virtual SigEst getSig() const {
    return SigEst();
  }

  /// Return the average non-diffractive impact parameter.
  double avNDB() const {
    return avNDb;
  }

  /// Calculate the Chi2 for the given cross section estimates.
  double Chi2(const SigEst & sigs, int npar) const;

  /// Use a simlified genetic algorithm to fit the parameters.
  virtual bool evolve();

  /// Set the parameters of this model.
  virtual void setParm(const vector<double> &) {}

  /// Return the current parameters and the minimum and maximum
  /// allowed values for the parameters of this model.
  virtual vector<double> getParm() const {
    return vector<double>();
  }
  virtual vector<double> minParm() const {
    return vector<double>();
  }
  virtual vector<double> maxParm() const {
    return vector<double>();
  }

private:

  /// The nucleon-nucleon cross sections targets for this model
  /// (tot, nd, dd, sdp, sdt, cd, el, bslope) and the required precision.
  vector<double> sigTarg, sigErr;

protected:

  /// The parameters stearing the fitting of internal parameters to
  /// the different nucleon-nucleon cross sections.
  int NInt, NGen, NPop;
  double sigFuzz;
  bool fitPrint;

  /// The estimated average impact parameter distance (in femtometer)
  /// for absorptive collisions.
  double avNDb;

  /// Info from the controlling HeavyIons object
  NucleusModel * projPtr;
  NucleusModel * targPtr;
  SigmaTotal * sigTotPtr;
  Settings * settingsPtr;
  Info * infoPtr;
  Rndm * rndPtr;

};


//==========================================================================

/// The most naive sub-collision model, asuming static nucleons and
/// the absorptive cross section equal to the total inelastic. No
///  fluctuations, meaning no diffraction.

class BlackSubCollisionModel: public SubCollisionModel {

public:

  /// The default constructor simply lists the nucleon-nucleon cross
  /// sections.
  BlackSubCollisionModel() {}

  /// Virtual destructor,
  virtual ~BlackSubCollisionModel() {}

  /// Take two vectors of Nucleons and an impact parameter vector and
  /// produce the corrsponding sub-collisions. Note that states of the
  /// nucleons may be changed.
  virtual multiset<SubCollision>
  getCollisions(vector<Nucleon> & proj, vector<Nucleon> & targ,
                const Vec4 & bvec, double & T);

};

//==========================================================================

/// A very simple sub-collision model, asuming static nucleons and
/// just assuring that the individual nucleon-nucleon cross sections
/// are preserved.

class NaiveSubCollisionModel: public SubCollisionModel {

public:

  /// The default constructor simply lists the nucleon-nucleon cross
  /// sections.
  NaiveSubCollisionModel() {}

  /// Virtual destructor,
  virtual ~NaiveSubCollisionModel() {}

  /// Take two vectors of Nucleons and an impact parameter vector and
  /// produce the corrsponding sub-collisions. Note that states of the
  /// nucleons may be changed.
  virtual multiset<SubCollision>
  getCollisions(vector<Nucleon> & proj, vector<Nucleon> & targ,
                const Vec4 & bvec, double & T);

};

//==========================================================================

/// A more complicated model where each nucleon has a fluctuating
/// "radius" according to a Strikman-inspired distribution.

class DoubleStrikman: public SubCollisionModel {

public:

  /// The default constructor simply lists the nucleon-nucleon cross
  /// sections.
  DoubleStrikman(int modein = 0)
    : r0(0.0), k0(4.0), sigd(75.0), alpha(0.5), opacityMode(modein) {}

  /// Virtual destructor,
  virtual ~DoubleStrikman() {}

  /// Take two vectors of Nucleons and an impact parameter vector and
  /// produce the corrsponding sub-collisions. Note that states of the
  /// nucleons may be changed.
  virtual multiset<SubCollision>
  getCollisions(vector<Nucleon> & proj, vector<Nucleon> & targ,
                const Vec4 & bvec, double & T);

  /// Generate a random number according to a Gamma-distribution.
  double gamma() const;

  /// The opacity of the collision at a given sigma.
  double opacity(double sig) const {
    // *** THINK *** maybe sig/sigd?
    sig /= sigd;
    if ( opacityMode == 1 ) sig = 1.0/sig;
    return sig > std::numeric_limits<double>::epsilon()?
      pow(-expm1(-1.0/sig), alpha): 1.0;
  }

  /// Return the elastic amplitude for a projectile and target state
  /// and the impact parameter between the corresponding nucleons.
  double Tpt(const Nucleon::State & p,
             const Nucleon::State & t, double b) const {
    double sig = M_PI*pow2(p[0] + t[0]);
    double grey = opacity(sig);
    return sig/grey > b*b*2.0*M_PI? grey: 0.0;
  }

  /// Calculate the cross sections for the given set of parameters.
  SigEst getSig() const;

  /// Set the parameters of this model.
  virtual void setParm(const vector<double> &);

  /// Return the current parameters and the minimum and maximum
  /// allowed values for the parameters of this model.
  virtual vector<double> getParm() const;
  virtual vector<double> minParm() const;
  virtual vector<double> maxParm() const;

  // Helper functions
  static void shuffle(double PND1, double PND2,
                      double & PW1, double & PW2);
  static void shuffel(double & PEL11, double P11,
                      double P12, double P21, double P22);
  static double PNW(double PWp, double PWt, double PND) {
    return ( 1.0 - PWp <= 0.0 || 1.0 - PWt <= 0.0 )?
      0.0: (1.0 - PWp)*(1.0 - PWt)/(1.0 - PND);
  }

protected:

  /// The average radius of the nucleon.
  double r0;

  /// The power in the Gamma distribution.
  double k0;

  /// Saturation scale of the nucleus.
  double sigd;

  /// Power of the saturation scale
  double alpha;

  /// Optional mode for opacity.
  int opacityMode;

};


//==========================================================================

/// A more complicated model where each nucleon has a fluctuating
/// "radius" according to a Strikman-inspired distribution.

class MultiRadial: public SubCollisionModel {

public:

  /// The default constructor simply lists the nucleon-nucleon cross
  /// sections.
  MultiRadial(int NrIn = 0)
    : Nr(max(1, NrIn)) {
    dR = T0 = c = phi = vector<double>(Nr, 0.0);
  }

  /// Virtual destructor,
  virtual ~MultiRadial() {}

  /// Take two vectors of Nucleons and an impact parameter vector and
  /// produce the corrsponding sub-collisions. Note that states of the
  /// nucleons may be changed.
  virtual multiset<SubCollision>
  getCollisions(vector<Nucleon> & proj, vector<Nucleon> & targ,
                const Vec4 & bvec, double & T);

  /// Return the elastic amplitude for a projectile and target state
  /// and the impact parameter between the corresponding nucleons.
  double Tpt(const Nucleon::State & p,
             const Nucleon::State & t, double b) const {
    return b < p[0] + t[0]? p[1]*t[1]: 0.0;
  }

  /// Calculate the cross sections for the given set of parameters.
  SigEst getSig() const;

  /// Set the parameters of this model.
  virtual void setParm(const vector<double> &);

  /// Return the current parameters and the minimum and maximum
  /// allowed values for the parameters of this model.
  virtual vector<double> getParm() const;
  virtual vector<double> minParm() const;
  virtual vector<double> maxParm() const;

protected:

  // Set the probabilities according to the angle parameters.
  void setProbs();

  /// Choose a radius.
  int choose() const;

  /// The number of radii.
  int Nr;


  /// The probability distribution.
  vector<double> c;

  /// The difference between radii.
  vector<double> dR;

  /// The opacity for different radii.
  vector<double> T0;

  /// The angles defining the probability distribution for the radii.
  vector<double> phi;

};


//==========================================================================

// Class for storing Events and Info objects.

class EventInfo {

public:

  /// Empty constructor.
  EventInfo(): ordering(-1.0), coll(0), ok(false) {}

  /// The Event and info objects.
  Event event;
  Info info;

  /// The ordering variable of this event.
  double ordering;
  bool operator<(const EventInfo & ei) const {
    return ordering < ei.ordering;
  }

  /// The associated SubCollision object.
  const SubCollision * coll;

  /// Is the event properly generated?
  bool ok;

  /// Which projectile and target nucleons are included and where are
  /// they placed?
  map<Nucleon *, pair<int,int> > projs, targs;

};

//==========================================================================

/// Class for collecting info about a Heavy Ion run and its produced
/// events.

class HIInfo {

public:

  friend class HeavyIons;
  friend class Angantyr;

  /// Constructor.
  HIInfo()
    : idProjSave(0), idTargSave(0), bSave(0.0), NSave(0), NAccSave(0),
      sigmaTotSave(0.0), sigmaNDSave(0.0), sigErr2TotSave(0.0),
      sigErr2NDSave(0.0), weightSave(0.0), weightSumSave(0.0),
      nCollSave(10, 0), nProjSave(10, 0), nTargSave(10, 0), nFailSave(0),
      subColsPtr(NULL) {}

  /// The impact-parameter distance in the current event.
  double b() const {
    return bSave;
  }

  /// The Monte Carlo integrated total cross section in the current run.
  double sigmaTot() const {
    return sigmaTotSave/millibarn;
  }

  /// The estimated statistical error on sigmaTot().
  double sigmaTotErr() const {
    return sqrt(sigErr2TotSave/max(1.0, double(NSave)))/millibarn;
  }

  /// The Monte Carlo integrated non-diffractive cross section in the
  /// current run.
  double sigmaND() const {
    return sigmaNDSave/millibarn;
  }

  /// The estimated statistical error on sigmaND().
  double sigmaNDErr() const {
    return sqrt(sigErr2NDSave/max(1.0, double(NSave)));
  }

  /// The number of attempted impact parameter points.
  long nAttempts() const {
    return NSave;
  }

  /// The number of produced events.
  long nAccepted() const {
    return NAccSave;
  }

  /// The total number of separate sub-collisions.
  int nCollTot() const { return nCollSave[0]; }

  /// The number of separate non-diffractive sub collisions in the
  /// current event.
  int nCollND() const { return nCollSave[1]; }

  /// The total number of non-diffractive sub collisions in the current event.
  int nCollNDTot() const { return nProjSave[1] + nTargSave[1] - nCollSave[1]; }

  /// The number of separate single diffractive projectile excitation
  /// sub collisions in the current event.
  int nCollSDP() const { return nCollSave[2]; }

  /// The number of separate single diffractive target excitation sub
  /// collisions in the current event.
  int nCollSDT() const { return nCollSave[3]; }

  /// The number of separate double diffractive sub collisions in the
  /// current event.
  int nCollDD() const { return nCollSave[4]; }

  /// The number of separate double diffractive sub collisions in the
  /// current event.
  int nCollCD() const { return nCollSave[5]; }

  /// The number of separate elastic sub collisions.
  int nCollEL() const { return nCollSave[6]; }

  /// The number of interacting projectile nucleons in the current
  /// event.
  int nPartProj() const { return nProjSave[0]; }

  /// The number of absorptively wounded projectile nucleons in the
  /// current event.
  int nAbsProj() const { return nProjSave[1]; }

  /// The number of diffractively wounded projectile nucleons in the
  /// current event.
  int nDiffProj() const { return nProjSave[2]; }

  /// The number of elastically scattered projectile nucleons in the
  /// current event.
  int nElProj() const { return nProjSave[3]; }

  /// The number of interacting projectile nucleons in the current
  /// event.
  int nPartTarg() const { return nTargSave[0]; }

  /// The number of absorptively wounded projectile nucleons in the
  /// current event.
  int nAbsTarg() const { return nTargSave[1]; }

  /// The number of diffractively wounded projectile nucleons in the
  /// current event.
  int nDiffTarg() const { return nTargSave[2]; }

  /// The number of elastically scattered projectile nucleons in the
  /// current event.
  int nElTarg() const { return nTargSave[3]; }

  /// The weight for this collision.
  double weight() const { return weightSave; }

  /// The sum of weights of the produced events.
  double weightSum() const { return weightSumSave; }

  /// The number of failed nuclon excitations in the current event.
  int nFail() const {
    return nFailSave;
  }

  /// Register a failed nucleon excitation.
  void failedExcitation() {
    ++nFailSave;
  }

 /// The number of separate non-diffractive collisions.
private:

  /// Register a tried impact parameter point giving the total elastic
  /// amplitude, the impact parameter and impact parameter generation weight.
  void addAttempt(double T, double bin, double bweight);

  /// Reweight event for whatever reason.
  void reweight(double w) {
    weightSave *= w;
  }

  // Select the primary process.
  void select(Info & in) {
    primInfo = in;
    primInfo.hiinfo = this;
  }

  // Accept an event and update statistics in info.
  void accept();

  /// Reject an attmpted event.
  void reject() {}

  /// Register a full sub collision.
  int addSubCollision(const SubCollision & c);

  /// Register a participating projectile/target nucleon.
  int addProjectileNucleon(const Nucleon & n);
  int addTargetNucleon(const Nucleon & n);


  /// Id of the colliding nuclei.
  int idProjSave, idTargSave;

  /// Impact parameter.
  double bSave;

  /// Cross section estimates.
  long NSave, NAccSave;
  double sigmaTotSave, sigmaNDSave, sigErr2TotSave, sigErr2NDSave;
  double weightSave, weightSumSave;

  /// Number of collisions and paricipants. See accessor functions for
  /// indices.
  vector<int> nCollSave, nProjSave, nTargSave;

  // Map of primary processes and the number of events and the sum of
  // weights.
  map<int,double> sumPrimW, sumPrimW2;
  map<int,int> NPrim;
  map<int,string> NamePrim;

  // The info object of the primary process.
  Info primInfo;

  // Number of failed nucleon excitations.
  int nFailSave;


public:
  // Access to subcollision to be extracted by the user.
  multiset<SubCollision>* subCollisionsPtr() { return subColsPtr; }

  void subCollisionsPtr(multiset<SubCollision> * sPtrIn) {
    subColsPtr = sPtrIn; }

private:

  // Full information about the Glauber calculation, consisting of
  // all subcollisions.
  multiset<SubCollision>* subColsPtr;

};

//==========================================================================

/// This is the heavy ion user hooks class which in the future may be
/// used inside a Pythia object to generate heavy ion collisons. For
/// now it is used outside Pythia and requires access to a number of
/// Pythia objects.

class HIUserHooks {

public:

  /// The default constructor is empty.
  HIUserHooks(): idProjSave(0), idTargSave(0) {}

  /// Virtual destructor.
  virtual ~HIUserHooks() {}

  /// Initialize this user hook.
  virtual void init(int idProjIn, int idTargIn) {
    idProjSave = idProjIn;
    idTargSave = idTargIn;
  }

  /// A user-supplied impact parameter generator.
  virtual bool hasImpactParameterGenerator() const { return false; }
  virtual ImpactParameterGenerator * impactParameterGenerator() const {
    return 0; }

  /// A suser-supplied NucleusModel for the projectile and target.
  virtual bool hasProjectileModel() const { return false; }
  virtual NucleusModel * projectileModel() const { return 0; }
  virtual bool hasTargetModel() const { return false; }
  virtual NucleusModel * targetModel() const { return 0; }

  /// A user-supplied SubCollisionModel for generating nucleon-nucleon
  /// subcollisions.
  virtual bool hasSubCollisionModel() { return false; }
  virtual SubCollisionModel * subCollisionModel() { return 0; }

  /// A user-supplied ordering of events in (inverse) hardness.
  virtual bool hasEventOrdering() const { return false; }
  virtual double eventOrdering(const Event &, const Info &) { return -1; }

  /// A user-supplied method for fixing up proton-neutron mismatch in
  /// generated beams.
  virtual bool canFixIsoSpin() const { return false; }
  virtual bool fixIsoSpin(EventInfo &) { return false; }

  /// A user-supplied method for shifting the event in impact parameter space.
  virtual bool canShiftEvent() const { return false; }
  virtual EventInfo & shiftEvent(EventInfo & ei) const { return ei; }

  /// A user-supplied method of adding a diffractive excitation event
  /// to another event, optionally connecting their colours.
  bool canAddNucleonExcitation() const { return false; }
  bool addNucleonExcitation(EventInfo &, EventInfo &, bool) const {
    return false; }

  /// A user supplied wrapper around the Pythia::forceHadronLevel()
  virtual bool canForceHadronLevel() const { return false; }
  virtual bool forceHadronLevel(Pythia &) { return false; }

  /// A user-supplied way of finding the remnants of an
  /// non-diffrcative pp collision (on the target side if tside is
  /// true) to be used to give momentum when adding.
  virtual bool canFindRecoilers() const { return false; }
  virtual vector<int>
  findRecoilers(const Event &, bool /* tside */, int /* beam */, int /* end */,
               const Vec4 & /* pdiff */, const Vec4 & /* pbeam */) const {
    return vector<int>();
  }

protected:

  /// Information set in the init() function.
  /// The PDG id of the projectile and target nuclei.
  int idProjSave, idTargSave;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HIUserHooks_H
