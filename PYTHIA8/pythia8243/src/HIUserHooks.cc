// HIUserHooks.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the HIUserHooks.h header) for
// the heavy ion classes classes, and some related global functions.

#include "Pythia8/Pythia.h"
#include "Pythia8/HIUserHooks.h"

namespace Pythia8 {


//==========================================================================

// NucleusModel base class.

//--------------------------------------------------------------------------

// Initialise base class, passing pointers to important objects.

void NucleusModel::initPtr(int idIn, Settings & settingsIn,
                           ParticleData & particleDataIn, Rndm & rndIn) {
  idSave = idIn;
  settingsPtr = &settingsIn;
  particleDataPtr = &particleDataIn;
  rndPtr = &rndIn;
  int decomp = abs(idSave);
  ISave = decomp%10;
  decomp /= 10;
  ASave = decomp%1000;
  decomp /= 1000;
  ZSave = decomp%1000;
  decomp /= 1000;
  LSave = decomp%10;
  decomp /= 10;

  if ( decomp != 10 ) {
    LSave = 0;
    ISave = 0;
    ASave = 0;
    ZSave = 0;
  }
}

//--------------------------------------------------------------------------

// Initialise in a subclass. Only a dummy in this base class.

bool NucleusModel::init() {
  return true;
}

//--------------------------------------------------------------------------

// Produce a proper particle corresponding to the nucleus handled by
// this NucleusModel.

Particle NucleusModel::produceIon(bool istarg) {
  double e = max(A(), 1)*(istarg? settingsPtr->parm("Beams:eB"):
                               settingsPtr->parm("Beams:eA"));
  double m = particleDataPtr->m0(id());
  Particle p(id(), -12);
  double pz = sqrt(max(e*e - m*m, 0.0));
  int daughter = 3;
  if ( istarg ) {
    pz = -pz;
    daughter = 4;
  }
  p.p(0.0, 0.0, pz, e);
  p.m(m);
  p.daughter1(daughter);
  return p;
}

//==========================================================================

// WoodsSaxonModel is a subclass of NucleusModel and implements a
// general Wood-Saxon distributed nucleus.

//--------------------------------------------------------------------------

// Place a nucleon inside a nucleus.

Vec4 WoodsSaxonModel::generateNucleon() const {

  while ( true ) {
    double r = R();
    double sel = rndPtr->flat()*(intlo + inthi0 + inthi1 + inthi2);
    if ( sel > intlo ) r -= a()*log(rndPtr->flat());
    if ( sel > intlo + inthi0 ) r -= a()*log(rndPtr->flat());
    if ( sel > intlo + inthi0 + inthi1 )  r -= a()*log(rndPtr->flat());
    if ( sel <= intlo ) {
      r = R()*pow(rndPtr->flat(), 1.0/3.0);
      if ( rndPtr->flat()*(1.0 + exp((r - R())/a())) > 1.0 ) continue;
    } else
      if ( rndPtr->flat()*(1.0 + exp((r - R())/a())) > exp((r - R())/a()) )
        continue;

    double costhe = 2.0*rndPtr->flat() - 1.0;
    double sinthe = sqrt(max(1.0 - costhe*costhe, 0.0));
    double phi = 2.0*M_PI*rndPtr->flat();

    return Vec4(r*sinthe*cos(phi), r*sinthe*sin(phi), r*costhe);

  }

  return Vec4();

}

//==========================================================================

// GLISSANDOModel is a subclass of NucleusModel and implements a
// general Wood-Saxon distributed nucleus with hard cores and specialy
// fitted parameters.

//--------------------------------------------------------------------------

// Initialize parameters.

bool GLISSANDOModel::init() {
  if ( A() == 0 ) return true;
  gaussHardCore = settingsPtr->flag("HeavyIon:gaussHardCore");
  if ( settingsPtr->isFlag("HI:hardCore") ) {
    if ( settingsPtr->flag("HI:hardCore") ) {
      RhSave = 0.9*femtometer;
      RSave = (1.1*pow(double(A()),1.0/3.0) -
               0.656*pow(double(A()),-1.0/3.0))*femtometer;
      aSave = 0.459*femtometer;
    } else {
      RSave = (1.12*pow(double(A()),1.0/3.0) -
               0.86*pow(double(A()),-1.0/3.0))*femtometer;
      aSave = 0.54*femtometer;
    }
    return WoodsSaxonModel::init();
  }
  if ( settingsPtr->flag("HeavyIon:WSHardCore") ) {
    RhSave = settingsPtr->parm("HeavyIon:WSRh");
    RSave = (1.1*pow(double(A()),1.0/3.0) -
             0.656*pow(double(A()),-1.0/3.0))*femtometer;
    aSave = 0.459*femtometer;
  } else {
    RSave = (1.12*pow(double(A()),1.0/3.0) -
             0.86*pow(double(A()),-1.0/3.0))*femtometer;
    aSave = 0.54*femtometer;
  }
  if ( settingsPtr->parm("HeavyIon:WSR") > 0.0 )
      RSave = settingsPtr->parm("HeavyIon:WSR");
  if ( settingsPtr->parm("HeavyIon:WSa") > 0.0 )
      aSave = settingsPtr->parm("HeavyIon:WSa");

  return WoodsSaxonModel::init();

}

//--------------------------------------------------------------------------

// Place a nucleon inside a nucleus.

vector<Nucleon> GLISSANDOModel::generate() const {
  int sign = id() > 0? 1: -1;
  int pid = sign*2212;
  int nid = sign*2112;
  vector<Nucleon> nucleons;
  if ( A() == 0 ) {
    nucleons.push_back(Nucleon(id(), 0, Vec4()));
    return nucleons;
  }
  if ( A() == 1 ) {
    if ( Z() == 1 ) nucleons.push_back(Nucleon(pid, 0, Vec4()));
    else  nucleons.push_back(Nucleon(nid, 0, Vec4()));
    return nucleons;
  }

  Vec4 cms;
  vector<Vec4> positions;
  while ( int(positions.size()) < A() ) {
    while ( true ) {
      Vec4 pos = generateNucleon();
      bool overlap = false;
      for ( int i = 0, N = positions.size(); i < N && !overlap; ++i )
        if ( (positions[i] - pos).pAbs() < (gaussHardCore ? RhGauss() : Rh()) )
          overlap = true;
      if ( overlap ) continue;
      positions.push_back(pos);
      cms += pos;
      break;
    }
  }

  cms /= A();
  nucleons.resize(A());
  int Np = Z();
  int Nn = A() - Z();
  for ( int i = 0, N= positions.size(); i < N; ++i ) {
    Vec4 pos(positions[i].px() - cms.px(),
                 positions[i].py() - cms.py());
    if ( int(rndPtr->flat()*(Np + Nn)) >= Np ) {
      --Nn;
      nucleons[i] = Nucleon(nid, i, pos);
    } else {
      --Np;
      nucleons[i] = Nucleon(pid, i, pos);
    }
  }

  return nucleons;

}

//==========================================================================

// ImpactParameterGenerator samples the impact parameter space.

//--------------------------------------------------------------------------

// Initialise base class, passing pointers to important objects.

void ImpactParameterGenerator::initPtr(SubCollisionModel & collIn,
                                       NucleusModel & projIn,
                                       NucleusModel & targIn,
                                       Settings & settingsIn,
                                       Rndm & rndIn) {
  collPtr = &collIn;
  projPtr = &projIn;
  targPtr = &targIn;
  settingsPtr = &settingsIn;
  rndPtr = &rndIn;

}

//--------------------------------------------------------------------------

// Initialise base class, bay be overridden by subclasses.

bool ImpactParameterGenerator::init() {
  if ( settingsPtr->isParm("HI:bWidth") )
    widthSave = settingsPtr->parm("HI:bWidth")*femtometer;
  else
    widthSave = settingsPtr->parm("HeavyIon:bWidth")*femtometer;

  if ( widthSave <= 0.0 ) {
    double Rp = sqrt(collPtr->sigTot()/M_PI)/2.0;
    double RA = max(Rp, projPtr->R());
    double RB = max(Rp, targPtr->R());
    widthSave = RA + RB + 2.0*Rp;
    cout << " HeavyIon Info: Initializing impact parameter generator "
         << "with width " << widthSave/femtometer << " fm." << endl;
  }

  return true;
}

//--------------------------------------------------------------------------

// Generate an impact parameter according to a gaussian distribution.

Vec4 ImpactParameterGenerator::generate(double & weight) const {
  double b = sqrt(-2.0*log(rndPtr->flat()))*width();
  double phi = 2.0*M_PI*rndPtr->flat();
  weight = 2.0*M_PI*width()*width()*exp(0.5*b*b/(width()*width()));
  return Vec4(b*sin(phi), b*cos(phi), 0.0, 0.0);
}

//==========================================================================

// The SubCollisionModel base class for modeling the collision between
// two nucleons to tell which type of collision has occurred. The
// model may manipulate the corresponing state of the nucleons.

//--------------------------------------------------------------------------

// Initialize the base class. Subclasses should consider calling this
// in overriding functions.

bool SubCollisionModel::init() {
  sigTarg[0] = sigTotPtr->sigmaTot()*millibarn;
  sigTarg[1] = sigTotPtr->sigmaND()*millibarn;
  sigTarg[2] = sigTotPtr->sigmaXX()*millibarn;
  sigTarg[3] = sigTotPtr->sigmaAX()*millibarn + sigTarg[1] + sigTarg[2];
  sigTarg[4] = sigTotPtr->sigmaXB()*millibarn + sigTarg[1] + sigTarg[2];
  sigTarg[5] = sigTotPtr->sigmaAXB()*millibarn;
  sigTarg[6] = sigTotPtr->sigmaEl()*millibarn;
  sigTarg[7] = sigTotPtr->bSlopeEl();
  NInt = settingsPtr->mode("HeavyIon:SigFitNInt");
  NGen = settingsPtr->mode("HeavyIon:SigFitNGen");
  NPop = settingsPtr->mode("HeavyIon:SigFitNPop");
  sigErr = settingsPtr->pvec("HeavyIon:SigFitErr");
  sigFuzz = settingsPtr->parm("HeavyIon:SigFitFuzz");
  fitPrint = settingsPtr->flag("HeavyIon:SigFitPrint");
  // preliminarily set average non-diffractive impact parameter as if
  // black disk.
  avNDb = 2.0*sqrt(sigTarg[1]/M_PI)*
    settingsPtr->parm("Angantyr:impactFudge")/3.0;
  return evolve();
}

//--------------------------------------------------------------------------

// Calculate the Chi^2 for the cross section that model in a subclass
// tries to mdoel.

double SubCollisionModel::Chi2(const SigEst & se, int npar) const {

  double chi2 = 0.0;
  int nval = 0;
  for ( int i = 0, Nval = se.sig.size(); i < Nval; ++i ) {
    if ( sigErr[i] == 0.0 ) continue;
    ++nval;
    chi2 += pow2(se.sig[i] - sigTarg[i])/
      (se.dsig2[i] + pow2(sigTarg[i]*sigErr[i]));
  }
  return chi2/double(max(nval - npar, 1));
}


//--------------------------------------------------------------------------

// Anonymous helper function to print out stuff.

namespace {

void printTarget(string name, double sig, double sigerr,
                 string unit = "mb    ") {
  cout << fixed << setprecision(2);
  cout << " |" << setw(25) << name << ": " << setw(8) << sig << " " << unit;
  if ( sigerr > 0.0 )
    cout <<"  (+- " << setw(2) << int(100.0*sigerr)
         << "%)                 | \n";
  else
    cout << "  not used                 | \n";
}

void printFit(string name, double fit, double sig, double sigerr,
                 string unit = "mb    ") {
  cout << " |" << setw(25) << name << ": "
       << setw(8) << fit
       << (sigerr > 0.0? " *(": "  (")
       << setw(6) << sig
       << ") " << unit << "                 | " << endl;
}

}
//--------------------------------------------------------------------------

// A simple genetic algorithm for fitting the parameters in a subclass
// to reproduce desired cross sections.

bool SubCollisionModel::evolve() {
  typedef vector<double> Parms;
  Parms minp = minParm();
  Parms maxp = maxParm();
  int dim = minp.size();
  if ( dim == 0 ) return true;

  if ( fitPrint ) {
    cout << " *------ HeavyIon fitting of SubCollisionModel to "
         << "cross sections ------* " << endl;
    printTarget("Total", sigTarg[0]/millibarn, sigErr[0]);
    printTarget("non-diffractive", sigTarg[1]/millibarn, sigErr[1]);
    printTarget("XX diffractive", sigTarg[2]/millibarn, sigErr[2]);
    printTarget("wounded target (B)", sigTarg[3]/millibarn, sigErr[3]);
    printTarget("wounded projectile (A)", sigTarg[4]/millibarn, sigErr[4]);
    printTarget("AXB diffractive", sigTarg[5]/millibarn, sigErr[5]);
    printTarget("elastic", sigTarg[6]/millibarn, sigErr[6]);
    printTarget("elastic b-slope", sigTarg[7], sigErr[7], "GeV^-2");
  }
  // We're going to use a home-made genetic algorithm. We start by
  // creating a population of random parameter points.
  vector<Parms> pop(NPop, Parms(dim));
  vector<double> def = settingsPtr->pvec("HeavyIon:SigFitDefPar");
  if ( settingsPtr->isPVec("HI:SigFitDefPar") )
    def = settingsPtr->pvec("HI:SigFitDefPar");
  for ( int j = 0; j < dim; ++j )
    pop[0][j] = max(minp[j], min(def[j], maxp[j]));
  for ( int i = 1; i < NPop; ++i )
    for ( int j = 0; j < dim; ++j )
      pop[i][j] = minp[j] + rndPtr->flat()*(maxp[j] - minp[j]);

  // Now we evolve our population for a number of generations.
  for ( int igen = 0; igen < NGen; ++igen ) {
    // Calculate Chi2 for each parameter set and order them.
    multimap<double, Parms> chi2map;
    double chi2max = 0.0;
    double chi2sum = 0.0;
    for ( int i = 0; i < NPop; ++i ) {
      setParm(pop[i]);
      double chi2 = Chi2(getSig(), dim);
      chi2map.insert(make_pair(chi2, pop[i]));
      chi2max = max(chi2max, chi2);
      chi2sum += chi2;
    }

    // Keep the best one, and move the other closer to a better one or
    // kill them if they are too bad.
    multimap<double, Parms>::iterator it  = chi2map.begin();
    pop[0] = it->second;
    if ( fitPrint ) {
      if ( igen == 0 )
        cout << " |                                      "
             << "                               | \n"
             << " |   Using a genetic algorithm          "
             << "                               | \n"
             << " |               Generation    best Chi2/Ndf     "
             << "                      | \n";
      cout << " |" << setw(25) << igen << ":" << fixed << setprecision(2)
           << setw(9) << it->first
           << "                                  | " << endl;
    }

    for ( int i = 1; i < NPop; ++i ) {
      pop[i] = (++it)->second;
      if ( it->first > rndPtr->flat()*chi2max ) {
        // Kill this individual and create a new one.
        for ( int j = 0; j < dim; ++j )
          pop[i][j] = minp[j] + rndPtr->flat()*(maxp[j] - minp[j]);
      } else {
        // Pick one of the better parameter sets and move this closer.
        int ii = int(rndPtr->flat()*i);
        for ( int j = 0; j < dim; ++j ) {
          double d = pop[ii][j] - it->second[j];
          double pl = max(minp[j], min(it->second[j] - sigFuzz*d, maxp[j]));
          double pu = max(minp[j], min(it->second[j] +
                                       (1.0 + sigFuzz)*d, maxp[j]));
          pop[i][j] = pl + rndPtr->flat()*(pu - pl);
        }
      }
    }
  }
  setParm(pop[0]);
  SigEst se = getSig();
  double chi2 = Chi2(se, dim);
  avNDb = se.avNDb*settingsPtr->parm("Angantyr:impactFudge");
  if ( chi2 > 2.0 )
    infoPtr->errorMsg("HeavyIon Warning: Chi^2 in fitting sub-collision "
                      "model to cross sections was high.");
  if ( fitPrint ) {
    cout << fixed << setprecision(2);
    cout << " |                                      "
         << "                               | "
         << endl;
    cout << " |     Resulting cross sections (target value) "
         << "                        | "
         << endl;
    printFit("Total", se.sig[0]/millibarn,
             sigTarg[0]/millibarn, sigErr[0]);
    printFit("non-diffractive", se.sig[1]/millibarn,
             sigTarg[1]/millibarn, sigErr[1]);
    printFit("XX diffractive", se.sig[2]/millibarn,
             sigTarg[2]/millibarn, sigErr[2]);
    printFit("wounded target (B)", se.sig[3]/millibarn,
             sigTarg[3]/millibarn, sigErr[3]);
    printFit("wounded projectile (A)", se.sig[4]/millibarn,
             sigTarg[4]/millibarn, sigErr[4]);
    printFit("AXB diffractive", se.sig[5]/millibarn,
             sigTarg[5]/millibarn, sigErr[5]);
    printFit("elastic", se.sig[6]/millibarn,
             sigTarg[6]/millibarn, sigErr[6]);
    printFit("elastic b-slope", se.sig[7], sigTarg[7], sigErr[7], "GeV^-2");
    cout << " |                 Chi2/Ndf: "
         << setw(8) << chi2
         << "                                  | " << endl;
    cout << " |                                   "
         << "                                  | "
         << endl;
    cout << " |     Resulting parameters:         "
         << "                                  | "
         << endl;
    for ( int j = 0; j < dim; ++j )
      cout << " |" << setw(25) << j << ":" << setw(9) << pop[0][j]
           << "                                  | " << endl;
    cout << " |                                      "
         << "                               | "
         << endl;
    cout << " |     Resulting non-diffractive average impact parameter: "
         << "            | "
         << endl;
    cout << " |                      <b>:" << setw(9) << se.avNDb << " +- "
         << setw(4) << sqrt(se.davNDb2) << " fm                       | "
         << endl;

    cout << " |                                      "
         << "                               | "
         << endl;
    cout << " *--- End HeavyIon fitting of parameters in "
         << "nucleon collision model ---* "
         << endl << endl;
    if ( NGen > 0 ) {
      cout << "HeavyIon Info: To avoid refitting, use the following settings "
           << "for next run:\n  HeavyIon:SigFitNGen = 0\n  "
           << "HeavyIon:SigFitDefPar = "
           << pop[0][0];
      for ( int j = 1; j < dim; ++j ) cout << "," << pop[0][j];
      for ( int j = dim; j < 8; ++j ) cout << ",0.0";
      cout << endl;
    }
  }

  return true;

}

//--------------------------------------------------------------------------

// For a given impact parameter and vectors of nucleons in the
// colliding nuclei, return a list of possible nucleon-nucleon
// SubCollisions ordered in nucleon-nucleon impact parameter
// distance. Should be called in overriding function in subclasses to
// reset all Nucleon objects.

multiset<SubCollision> SubCollisionModel::
getCollisions(vector<Nucleon> & proj, vector<Nucleon> & targ,
              const Vec4 & bvec, double & T) {
  multiset<SubCollision> ret;
  T = 0.0;
  // Reset all states.
  for ( int i = 0, N = proj.size(); i < N; ++i ) {
    proj[i].reset();
    proj[i].bShift(bvec/2.0);
  }
  for ( int i = 0, N = targ.size(); i < N; ++i ) {
    targ[i].reset();
    targ[i].bShift(-bvec/2.0);
  }

  // Empty return.
  return ret;
}


//==========================================================================

// The BlackSubCollisionModel uses fixed size, black-disk
// nucleon-nucleon cross section, equal to the total inelastic pp cross
// section. Everything else is elastic -- Diffraction not included.

//--------------------------------------------------------------------------

multiset<SubCollision> BlackSubCollisionModel::
getCollisions(vector<Nucleon> & proj, vector<Nucleon> & targ,
              const Vec4 & bvec, double & T) {
  // Always call base class to reset nucleons and shift them into
  // position.
  multiset<SubCollision> ret =
    SubCollisionModel::getCollisions(proj, targ, bvec, T);

  T = 0.0;
  // Go through all pairs of nucleons
  for ( int ip = 0, Np = proj.size(); ip < Np; ++ip )
    for ( int it = 0, Nt = targ.size(); it < Nt; ++it ) {
      Nucleon & p = proj[ip];
      Nucleon & t = targ[it];
      double b = (p.bPos() - t.bPos()).pT();
      if ( b > sqrt(sigTot()/M_PI) ) continue;
      T = 0.5; // The naive cross section only gets the total xsec correct.
      if ( b < sqrt((sigTot() - sigEl())/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ABS));
      }
      else {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ELASTIC));
      }
    }

  return ret;
}

//==========================================================================

// The NaiveSubCollisionModel uses a fixed size, black-disk-like
// nucleon-nucleon cross section where. Central collisions will always
// be absorptive, less central will be doubly diffractive, more
// peripheral will be single diffractive and the most peripheral will
// be elastic.

//--------------------------------------------------------------------------

multiset<SubCollision> NaiveSubCollisionModel::
getCollisions(vector<Nucleon> & proj, vector<Nucleon> & targ,
              const Vec4 & bvec, double & T) {
  // Always call base class to reset nucleons and shift them into
  // position.
  multiset<SubCollision> ret =
    SubCollisionModel::getCollisions(proj, targ, bvec, T);

  T = 0.0;
  // Go through all pairs of nucleons
  for ( int ip = 0, Np = proj.size(); ip < Np; ++ip )
    for ( int it = 0, Nt = targ.size(); it < Nt; ++it ) {
      Nucleon & p = proj[ip];
      Nucleon & t = targ[it];
      double b = (p.bPos() - t.bPos()).pT();
      if ( b > sqrt(sigTot()/M_PI) ) continue;
      T = 0.5; // The naive cross section only gets the total xsec correct.
      if ( b < sqrt(sigND()/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ABS));
      }
      else if ( b < sqrt((sigND() + sigDDE())/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::DDE));
      }
      else if ( b < sqrt((sigND() + sigSDE() + sigDDE())/M_PI) ) {
         if ( sigSDEP() > rndPtr->flat()*sigSDE() ) {
          ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::SDEP));
        } else {
          ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::SDET));
        }
      }
      else if ( b < sqrt((sigND() + sigSDE() + sigDDE() + sigCDE())/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::CDE));
      }
      else {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ELASTIC));
      }
    }

  return ret;
}

//==========================================================================

// DoubleStrikman uses a fluctuating and semi-transparent disk for
// each Nucleon in a sub-collision resulting in a fluctuating
// interaction probability. To assess the fluctuation each Nucleon has
// two random states in each collision, one main state and one helper
// state to assess the frlutuations.

//--------------------------------------------------------------------------

// Anonymous helper functions to simplify calculating elastic
// amplitudes.

namespace {
inline double el(double s1, double s2, double u1, double u2) {
  return s1/u1 > s2/u2? s2*u1: s1*u2;
}

}

//--------------------------------------------------------------------------

// Numerically estimate the cross sections corresponding to the
// current parameter setting.

SubCollisionModel::SigEst DoubleStrikman::getSig() const {

  SigEst s;
  for ( int n = 0; n < NInt; ++n ) {
    double rp1 = gamma();
    double rp2 = gamma();
    double rt1 = gamma();
    double rt2 = gamma();
    double s11 = pow2(rp1 + rt1)*M_PI;
    double s12 = pow2(rp1 + rt2)*M_PI;
    double s21 = pow2(rp2 + rt1)*M_PI;
    double s22 = pow2(rp2 + rt2)*M_PI;

    double stot = (s11 + s12 + s21 + s22)/4.0;
    s.sig[0] += stot;
    s.dsig2[0] += pow2(stot);

    double u11 = opacity(s11)/2.0;
    double u12 = opacity(s12)/2.0;
    double u21 = opacity(s21)/2.0;
    double u22 = opacity(s22)/2.0;

    double avb = sqrt(2.0/M_PI)*(s11*sqrt(s11/(2.0*u11))*(1.0 - u11) +
                                 s12*sqrt(s12/(2.0*u12))*(1.0 - u12) +
                                 s21*sqrt(s21/(2.0*u21))*(1.0 - u21) +
                                 s22*sqrt(s22/(2.0*u22))*(1.0 - u22))/12.0;
    s.avNDb += avb;
    s.davNDb2 += pow2(avb);

    double snd = (s11 - s11*u11 + s12 - s12*u12 +
                  s21 - s21*u21 + s22 - s22*u22)/4.0;
    s.sig[1] += snd;
    s.dsig2[1] += pow2(snd);

    double sel = (el(s11, s22, u11, u22) + el(s12, s21, u12, u21))/2.0;
    s.sig[6] += sel;
    s.dsig2[6] += pow2(sel);

    double swt = stot - (el(s11, s12, u11, u12) + el(s21, s22, u21, u22))/2.0;
    double swp = stot - (el(s11, s21, u11, u21) + el(s12, s22, u12, u22))/2.0;
    s.sig[4] += swp;
    s.dsig2[4] += pow2(swp);
    s.sig[3] += swt;
    s.dsig2[3] += pow2(swt);

    s.sig[2] += swt + swp - snd  + sel - stot;
    s.dsig2[2] += pow2(swt + swp - snd  + sel - stot);

    s.sig[5] += s11;
    s.dsig2[5] += pow2(s11);

    s.sig[7] += pow2(s11)/u11;
    s.dsig2[7] += pow2(pow2(s11)/u11);



  }

  s.sig[0] /= double(NInt);
  s.dsig2[0] = (s.dsig2[0]/double(NInt) - pow2(s.sig[0]))/double(NInt);

  s.sig[1] /= double(NInt);
  s.dsig2[1] = (s.dsig2[1]/double(NInt) - pow2(s.sig[1]))/double(NInt);

  s.sig[2] /= double(NInt);
  s.dsig2[2] = (s.dsig2[2]/double(NInt) - pow2(s.sig[2]))/double(NInt);

  s.sig[3] /= double(NInt);
  s.dsig2[3] = (s.dsig2[3]/double(NInt) - pow2(s.sig[3]))/double(NInt);

  s.sig[4] /= double(NInt);
  s.dsig2[4] = (s.dsig2[4]/double(NInt) - pow2(s.sig[4]))/double(NInt);

  s.sig[6] /= double(NInt);
  s.dsig2[6] = (s.dsig2[6]/double(NInt) - pow2(s.sig[6]))/double(NInt);

  s.sig[5] /= double(NInt);
  s.dsig2[5] /= double(NInt);

  s.sig[7] /= double(NInt);
  s.dsig2[7] /= double(NInt);
  double bS = (s.sig[7]/s.sig[5])/(16.0*M_PI*pow2(0.19732697));
  double b2S = pow2(bS)*(s.dsig2[7]/pow2(s.sig[7]) - 1.0 +
                        s.dsig2[5]/pow2(s.sig[5]) - 1.0)/double(NInt);
  s.sig[5] = 0.0;
  s.dsig2[5] = 0.0;
  s.sig[7] = bS;
  s.dsig2[7] = b2S;

  s.avNDb /= double(NInt);
  s.davNDb2 = (s.davNDb2/double(NInt) - pow2(s.avNDb))/double(NInt);
  s.avNDb /= s.sig[1];
  s.davNDb2 /= pow2(s.sig[1]);

  return s;

}

//--------------------------------------------------------------------------

// Access funtions to parameters in the DoubleStrikman model.

void DoubleStrikman::setParm(const vector<double> & p) {
  if ( p.size() > 0 ) sigd = p[0];
  if ( p.size() > 1 ) k0 = p[1];
  if ( p.size() > 2 ) alpha = p[2];
  r0 = sqrt(sigTot()/(M_PI*(2.0*k0 + 4.0*k0*k0)));
}

vector<double> DoubleStrikman::getParm() const {
  vector<double> ret(3);
  ret[0] = sigd;
  ret[1] = k0;
  ret[2] = alpha;
  return ret;
}

vector<double> DoubleStrikman::minParm() const {
  vector<double> ret(3);
  ret[0] = 1.0;
  ret[1] = 0.01;
  ret[2] = 0.0;
  return ret;
}

vector<double> DoubleStrikman::maxParm() const {
  vector<double> ret(3);
  ret[0] = 20.0;
  ret[1] = 20.0;
  ret[2] = 2.0;
  return ret;
}

//--------------------------------------------------------------------------

// Generate the radii in DoubleStrikman according to a gamma
// distribution.

double DoubleStrikman::gamma() const {
  static const double e = exp(1);
  int k = int(k0);
  double del = k0 - k;
  double x = 0.0;
  for ( int i = 0; i < k; ++i ) x += -log(rndPtr->flat());

  if ( del == 0.0 ) return x*r0;

  while ( true ) {
    double U = rndPtr->flat();
    double V = rndPtr->flat();
    double W = rndPtr->flat();

    double xi = 0.0;
    if ( U <= e/(e+del) ) {
      xi = pow(V, 1.0/del);
      if ( W <= exp(-xi) ) return r0*(x + xi);
    } else {
      xi = 1.0 - log(V);
      if ( W <= pow(xi, del - 1.0) ) return r0*(x + xi);
    }
  }
  return 0.0;
}

//--------------------------------------------------------------------------

// Helper functions to get the correct average elastic and wounded
// cross sections.

void DoubleStrikman::shuffle(double PND1, double PND2,
                             double & PW1, double & PW2) {
  if ( PND1 > PW1 ) {
    PW2 += PW1 - PND1;
    PW1 = PND1;
    return;
  }
  if ( PND2 > PW2 ) {
    PW1 += PW2 - PND2;
    PW2 = PND2;
    return;
  }
}

void DoubleStrikman::shuffel(double & PEL11, double P11,
                             double P12, double P21, double P22) {
  double PEL12 = PEL11, PEL21 = PEL11, PEL22 = PEL11;
  map<double, double *> ord;
  ord[P11] = &PEL11;
  ord[P12] = &PEL12;
  ord[P21] = &PEL21;
  ord[P22] = &PEL22;
  map<double, double *>::iterator next = ord.begin();
  map<double, double *>::iterator prev = next++;
  while ( next != ord.end() ) {
    if ( *prev->second > prev->first ) {
      *next->second += *prev->second - prev->first;
      *prev->second = prev->first;
    }
    prev = next++;
  }
}

//--------------------------------------------------------------------------

// Main function returning the possible sub-collisions.

multiset<SubCollision> DoubleStrikman::
getCollisions(vector<Nucleon> & proj, vector<Nucleon> & targ,
              const Vec4 & bvec, double & T) {
  // Always call base class to reset nucleons and shift them into
  // position.
  multiset<SubCollision> ret =
    SubCollisionModel::getCollisions(proj, targ, bvec, T);

  /// Assign two states to each nucleon
  for ( int ip = 0, Np = proj.size(); ip < Np; ++ip ) {
    proj[ip].state(Nucleon::State(1, gamma()));
    proj[ip].addAltState(Nucleon::State(1, gamma()));
  }
  for ( int it = 0, Nt = targ.size(); it < Nt; ++it ) {
    targ[it].state(Nucleon::State(1, gamma()));
    targ[it].addAltState(Nucleon::State(1, gamma()));
  }

  // The factorising S-matrix.
  double S = 1.0;

  // Go through all pairs of nucleons
  for ( int ip = 0, Np = proj.size(); ip < Np; ++ip )
    for ( int it = 0, Nt = targ.size(); it < Nt; ++it ) {
      Nucleon & p = proj[ip];
      Nucleon & t = targ[it];
      double b = (p.bPos() - t.bPos()).pT();

      double T11 = Tpt(p.state(), t.state(), b);
      double T12 = Tpt(p.state(), t.altState(), b);
      double T21 = Tpt(p.altState(), t.state(), b);
      double T22 = Tpt(p.altState(), t.altState(), b);
      double S11 = 1.0 - T11;
      double S12 = 1.0 - T12;
      double S21 = 1.0 - T21;
      double S22 = 1.0 - T22;
      S *= S11;
      double PND11 = 1.0 - pow2(S11);
      // First and most important, check if this is an absorptive
      // scattering.
      if ( PND11 > rndPtr->flat() ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ABS));
        continue;
      }

      // Now set up calculation for probability of diffractively
      // wounded nucleons.
      double PND12 = 1.0 - pow2(S12);
      double PND21 = 1.0 - pow2(S21);
      double PWp11 = 1.0 - S11*S21;
      double PWp21 = 1.0 - S11*S21;
      shuffle(PND11, PND21, PWp11, PWp21);
      double PWt11 = 1.0 - S11*S12;
      double PWt12 = 1.0 - S11*S12;
      shuffle(PND11, PND12, PWt11, PWt12);

      bool wt = ( PWt11 - PND11 > (1.0 - PND11)*rndPtr->flat() );
      bool wp = ( PWp11 - PND11 > (1.0 - PND11)*rndPtr->flat() );
      if ( wt && wp ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::DDE));
        continue;
      }
      if ( wt ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::SDET));
        continue;
      }
      if ( wp ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::SDEP));
        continue;
      }

      // Finally set up calculation for elastic scattering. This can
      // never be exact, but let's do as well as we can.

      double PND22 = 1.0 - pow2(S22);
      double PWp12 = 1.0 - S12*S22;
      double PWp22 = 1.0 - S12*S22;
      shuffle(PND12, PND22, PWp12, PWp22);
      double PWt21 = 1.0 - S21*S22;
      double PWt22 = 1.0 - S21*S22;
      shuffle(PND21, PND22, PWt21, PWt22);

      double PNW11 = PNW(PWp11, PWt11, PND11);
      double PNW12 = PNW(PWp12, PWt12, PND12);
      double PNW21 = PNW(PWp21, PWt21, PND21);
      double PNW22 = PNW(PWp22, PWt22, PND22);

      double PEL = (T12*T21 + T11*T22)/2.0;
      shuffel(PEL, PNW11, PNW12, PNW21, PNW22);
      if ( PEL > PNW11*rndPtr->flat() ) {
        if ( sigCDE() > rndPtr->flat()*(sigCDE() + sigEl()) )
          ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::CDE));
        else
          ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ELASTIC));
      }
    }


  T = 1.0 - S;

  return ret;
}

//==========================================================================

// MultiRadial uses a number of different disk sizes with different
// opacities. Like a discrete version of DoubleStrikman.

//--------------------------------------------------------------------------

// Numerically estimate the cross sections corresponding to the
// current parameter setting.

SubCollisionModel::SigEst MultiRadial::getSig() const {

  SigEst s;

  double sTpt = 0.0;
  double sT2pt = 0.0;
  //  double sTpt2 = 0.0;
  double sTp2t = 0.0;
  //  double sTt2p = 0.0;
  double Rp = 0.0;
  for ( int ip = 0; ip < Nr; ++ip ) {
    Rp += dR[ip];
    double Rt = 0.0;
    for ( int it = 0; it < Nr; ++it ) {
      Rt += dR[it];
      sTpt += c[ip]*T0[ip]*c[it]*T0[it]*pow2(Rp + Rt)*sigTot();
      sT2pt += c[ip]*pow2(T0[ip])*c[it]*pow2(T0[it])*pow2(Rp + Rt)*sigTot();
      double rp = 0.0;
      for ( int jp = 0; jp < Nr; ++jp ) {
        rp += dR[jp];
        double rt = 0.0;
        for ( int jt = 0; jt < Nr; ++jt ) {
          rt += dR[jt];
          double fac = T0[ip]*T0[jp]*T0[it]*T0[jt]*pow2(min(Rp + Rt, rp + rt))
            * sigTot();
          if ( ip == jp ) sTp2t += c[ip]*c[it]*c[jt]*fac;
        }
      }
    }

  }

  s.sig[0] /= double(NInt);
  s.dsig2[0] = (s.dsig2[0]/double(NInt) - pow2(s.sig[0]))/double(NInt);

  s.sig[1] /= double(NInt);
  s.dsig2[1] = (s.dsig2[1]/double(NInt) - pow2(s.sig[1]))/double(NInt);

  s.sig[2] /= double(NInt);
  s.dsig2[2] = (s.dsig2[2]/double(NInt) - pow2(s.sig[2]))/double(NInt);

  s.sig[3] /= double(NInt);
  s.dsig2[3] = (s.dsig2[3]/double(NInt) - pow2(s.sig[3]))/double(NInt);

  s.sig[4] /= double(NInt);
  s.dsig2[4] = (s.dsig2[4]/double(NInt) - pow2(s.sig[4]))/double(NInt);

  s.sig[6] /= double(NInt);
  s.dsig2[6] = (s.dsig2[6]/double(NInt) - pow2(s.sig[6]))/double(NInt);

  s.sig[5] /= double(NInt);
  s.dsig2[5] /= double(NInt);

  s.sig[7] /= double(NInt);
  s.dsig2[7] /= double(NInt);
  double bS = (s.sig[7]/s.sig[5])/(16.0*M_PI*pow2(0.19732697));
  double b2S = pow2(bS)*(s.dsig2[7]/pow2(s.sig[7]) - 1.0 +
                        s.dsig2[5]/pow2(s.sig[5]) - 1.0)/double(NInt);
  s.sig[5] = 0.0;
  s.dsig2[5] = 0.0;
  s.sig[7] = bS;
  s.dsig2[7] = b2S;

  return s;

}

//--------------------------------------------------------------------------

// Access funtions to parameters in the MultiRadial model.

void MultiRadial::setParm(const vector<double> & p) {
  unsigned int ip = 0;
  for ( int i = 0; i < Nr; ++i ) {
    if ( ip < p.size() ) dR[i] = p[ip++];
    if ( ip < p.size() ) T0[i] = p[ip++];
    if ( ip < p.size() ) phi[i] = p[ip++];
  }
}

vector<double> MultiRadial::getParm() const {
  vector<double> ret;
  for ( int i = 0; i < Nr; ++i ) {
    ret.push_back(dR[i]);
    ret.push_back(T0[i]);
    if ( i < Nr -1 ) ret.push_back(phi[i]);
  }
  return ret;
}

vector<double> MultiRadial::minParm() const {
  return vector<double>(Nr*Nr*(Nr - 1), 0.0);
}

vector<double> MultiRadial::maxParm() const {
  return vector<double>(Nr*Nr*(Nr - 1), 1.0);
}

void MultiRadial::setProbs() {
  double rProj = 1.0;
  for ( int i = 0; i < Nr - 1; ++i ) {
    c[i] = rProj*cos(phi[i]*M_PI/2.0);
    rProj *= sin(phi[i]*M_PI/2.0);
  }
  c[Nr - 1] = rProj;
}

int MultiRadial::choose() const {
  double sum = 0.0;
  double sel = rndPtr->flat();
  for ( int i = 0; i < Nr - 1; ++i )
    if ( sel < ( sum += c[i] ) ) return i;
  return Nr - 1;
}




//--------------------------------------------------------------------------

// Main function returning the possible sub-collisions.

multiset<SubCollision> MultiRadial::
getCollisions(vector<Nucleon> & proj, vector<Nucleon> & targ,
              const Vec4 & bvec, double & T) {
  // Always call base class to reset nucleons and shift them into
  // position.
  multiset<SubCollision> ret =
    SubCollisionModel::getCollisions(proj, targ, bvec, T);

  return ret;
}

//==========================================================================

// HIInfo functions to collect statistics in an event and in a run.

//--------------------------------------------------------------------------

// Collect statistics for each SubCollision in an event.

int HIInfo::addSubCollision(const SubCollision & c) {
  ++nCollSave[0];
  switch ( c.type ) {
  case SubCollision::ABS:
    return ++nCollSave[1];
  case SubCollision::SDEP:
    return ++nCollSave[2];
  case SubCollision::SDET:
    return ++nCollSave[3];
  case SubCollision::DDE:
    return ++nCollSave[4];
  case SubCollision::CDE:
    return ++nCollSave[5];
  case SubCollision::ELASTIC:
    return ++nCollSave[6];
  default:
    return 0;
  }
}

//--------------------------------------------------------------------------

// Collect statistics for each participating Nucleon in an event.

int HIInfo::addProjectileNucleon(const Nucleon & n) {
  ++nProjSave[0];
  switch ( n.status() ) {
  case Nucleon::ABS:
    return ++nProjSave[1];
  case Nucleon::DIFF:
    return ++nProjSave[2];
  case Nucleon::ELASTIC:
    return ++nProjSave[3];
  default:
    return 0;
  }
}

int HIInfo::addTargetNucleon(const Nucleon & n) {
  ++nTargSave[0];
  switch ( n.status() ) {
  case Nucleon::ABS:
    return ++nTargSave[1];
  case Nucleon::DIFF:
    return ++nTargSave[2];
  case Nucleon::ELASTIC:
    return ++nTargSave[3];
  default:
    return 0;
  }
}

//--------------------------------------------------------------------------

// Collect statistics for attemted and accepted impact paramet point
// in an event.

void HIInfo::addAttempt(double T, double bin, double bweight) {
  bSave = bin;
  nCollSave = nProjSave = nTargSave = vector<int>(10, 0);
  nFailSave = 0;
  weightSave = bweight;
  weightSumSave += weightSave;
  ++NSave;
  double w = 2.0*T*bweight;
  double delta = w - sigmaTotSave;
  sigmaTotSave += delta/double(NSave);
  sigErr2TotSave += (delta*(w - sigmaTotSave) - sigErr2TotSave)/double(NSave);
  w = (2*T - T*T)*bweight;
  delta = w - sigmaNDSave;
  sigmaNDSave += delta/double(NSave);
  sigErr2NDSave += (delta*(w - sigmaNDSave) - sigErr2NDSave)/double(NSave);
}


void HIInfo::accept() {
  int pc = primInfo.code();
  weightSumSave += weightSave;
  ++NAccSave;
  sumPrimW[pc] += weightSave;
  sumPrimW2[pc] += weightSave*weightSave;
  ++NPrim[pc];
  NamePrim[pc] = primInfo.nameProc(pc);
}

//==========================================================================

// The Nucleon class describing a Nucleon in a Nucleus.

//--------------------------------------------------------------------------

// Print out information about a Nuclean. To be called from inside a
// debugger.

void Nucleon::debug() {
  cout << "Nucleon id: " << id() << endl;
  cout << "index:      " << index() << endl;
  cout << "b(rel):     " << nPos().px() << " " << nPos().py() << endl;
  cout << "b(abs):     " << bPos().px() << " " << bPos().py() << endl;
  cout << "status:     " << status() << (done()? " done": "     ") << endl;
  cout << "state:      ";
  for ( int i = 0, N = state().size(); i < N; ++i )
    cout << state()[i] << " ";
  cout << endl;
  for ( int j = 0, M = altStatesSave.size(); j < M; ++j ) {
    cout << "state " << j+1 << ":    ";
    for ( int i = 0, N = altStatesSave[j].size(); i < N; ++i )
      cout << altStatesSave[j][i] << " ";
    cout << endl;
  }
}

//==========================================================================

} // end namespace Pythia8
