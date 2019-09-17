// HelicityMatrixElements.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for a number of physics classes used in tau decays.

#ifndef Pythia8_HelicityMatrixElements_H
#define Pythia8_HelicityMatrixElements_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/HelicityBasics.h"
#include "Pythia8/PythiaComplex.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/StandardModel.h"

namespace Pythia8 {

//==========================================================================

// The helicity matrix element class.

class HelicityMatrixElement {

public:

  // Constructor and destructor.
  HelicityMatrixElement() : DECAYWEIGHTMAX(), particleDataPtr(),
    couplingsPtr(), settingsPtr() {};
  virtual ~HelicityMatrixElement() {};

  // Initialize the physics matrices and pointers.
  virtual void initPointers(ParticleData*, Couplings*, Settings* = 0);

  // Initialize the channel.
  virtual HelicityMatrixElement* initChannel(vector<HelicityParticle>&);

  // Calculate the matrix element weight for a decay.
  virtual double decayWeight(vector<HelicityParticle>&);

  // Calculate the maximum matrix element decay weight.
  virtual double decayWeightMax(vector<HelicityParticle>&)
    {return DECAYWEIGHTMAX;}

  // Calculate the helicity matrix element.
  virtual complex calculateME(vector<int>){return complex(0,0);}

  // Calculate the decay matrix for a particle.
  virtual void calculateD(vector<HelicityParticle>&);

  // Calculate the density matrix for a particle.
  virtual void calculateRho(unsigned int, vector<HelicityParticle>&);

  // Set a fermion line.
  void setFermionLine(int, HelicityParticle&, HelicityParticle&);

  // Calculate Breit-Wigner's with running widths and fixed.
  virtual complex  breitWigner(double s, double M, double G);
  virtual complex sBreitWigner(double m0, double m1, double s,
    double M, double G);
  virtual complex pBreitWigner(double m0, double m1, double s,
    double M, double G);
  virtual complex dBreitWigner(double m0, double m1, double s,
    double M, double G);

protected:

  // Maximum decay weight. WARNING: hardcoded constant.
  double DECAYWEIGHTMAX;

  // Physics matrices.
  vector< GammaMatrix > gamma;

  // Particle map vector.
  vector< int > pMap;

  // Particle ID vector.
  vector< int > pID;

  // Particle mass vector.
  vector< double > pM;

  // Wave functions.
  vector< vector< Wave4 > > u;

  // Initialize the constants for the matrix element (called by initChannel).
  virtual void initConstants() {};

  // Initialize the wave functions (called by decayWeight and calculateRho/D).
  virtual void initWaves(vector<HelicityParticle>&) {};

  // Pointer to particle data.
  ParticleData* particleDataPtr;

  // Pointer to Standard Model constants.
  Couplings*    couplingsPtr;

  // Pointer to Settings.
  Settings*     settingsPtr;

private:

  // Recursive sub-method to calculate the density matrix for a particle.
  void calculateRho(unsigned int, vector<HelicityParticle>&,
    vector<int>&, vector<int>&, unsigned int);

  // Recursive sub-method to calculate the decay matrix for a particle.
  void calculateD(vector<HelicityParticle>&, vector<int>&, vector<int>&,
    unsigned int);

  // Recursive sub-method to calculate the matrix element weight for a decay.
  void decayWeight(vector<HelicityParticle>&, vector<int>&, vector<int>&,
    complex&, unsigned int);

  // Calculate the product of the decay matrices for a hard process.
  complex calculateProductD(unsigned int, unsigned int,
    vector<HelicityParticle>&, vector<int>&, vector<int>&);

  // Calculate the product of the decay matrices for a decay process.
  complex calculateProductD(vector<HelicityParticle>&,
    vector<int>&, vector<int>&);

};

//==========================================================================

// Helicity matrix element for the hard process of two fermions -> W/W' ->
// two fermions.

class HMETwoFermions2W2TwoFermions : public HelicityMatrixElement {

public:

  HMETwoFermions2W2TwoFermions() : p0CA(), p2CA(), p0CV(), p2CV() {}

  void initConstants();

  void initWaves(vector<HelicityParticle>&);

  complex calculateME(vector<int>);

private:

  // Vector and axial couplings.
  double p0CA, p2CA, p0CV, p2CV;

};

//==========================================================================

// Helicity matrix element for the hard process of two fermions ->
// photon/Z/Z' -> two fermions.

class HMETwoFermions2GammaZ2TwoFermions : public HelicityMatrixElement {

public:

  HMETwoFermions2GammaZ2TwoFermions() : p0CAZ(), p2CAZ(), p0CVZ(), p2CVZ(),
    p0CAZp(), p2CAZp(), p0CVZp(), p2CVZp(), cos2W(), sin2W(), zG(), zM(),
    zpG(), zpM(), s(), p0Q(), p2Q(), zaxis(), includeGamma(), includeZ(),
    includeZp() {}

  void initConstants();

  void initWaves(vector<HelicityParticle>&);

  complex calculateME(vector<int>);

private:

  // Return gamma element for the helicity matrix element.
  complex calculateGammaME(vector<int>);

  // Return Z/Z' element for helicity matrix element.
  complex calculateZME(vector<int>, double, double, double, double,
    double, double);

  // Return the Z' vector or axial coupling for a fermion.
  double zpCoupling(int id, string type);

  // Vector and axial couplings.
  double p0CAZ, p2CAZ, p0CVZ, p2CVZ, p0CAZp, p2CAZp, p0CVZp, p2CVZp;

  // Weinberg angle, Z/Z' width (on-shell), Z/Z' mass (on-shell), and s.
  double cos2W, sin2W, zG, zM, zpG, zpM, s;

  // Fermion line charge.
  double p0Q, p2Q;

  // Bool whether the incoming fermions are oriented with the z-axis.
  bool zaxis, includeGamma, includeZ, includeZp;

};

//==========================================================================

// Helicity matrix element for the hard process of X -> two fermions.

class HMEX2TwoFermions : public HelicityMatrixElement {

public:

  void initWaves(vector<HelicityParticle>&);

};

//==========================================================================

// Helicity matrix element for the hard process of W/W' -> two fermions.

class HMEW2TwoFermions : public HMEX2TwoFermions {

public:

  HMEW2TwoFermions() : p2CA(), p2CV() {}

  void initConstants();

  complex calculateME(vector<int>);

private:

  // Vector and axial couplings.
  double p2CA, p2CV;

};

//==========================================================================

// Helicity matrix element for the hard process of photon -> two fermions.

class HMEGamma2TwoFermions : public HMEX2TwoFermions {

public:

  complex calculateME(vector<int>);

};

//==========================================================================

// Helicity matrix element for the hard process of Z/Z' -> two fermions.

class HMEZ2TwoFermions : public HMEX2TwoFermions {

public:

  HMEZ2TwoFermions() : p2CA(), p2CV() {}

  void initConstants();

  complex calculateME(vector<int>);

private:

  // Return the Z' vector or axial coupling for a fermion.
  double zpCoupling(int id, string type);

  // Vector and axial couplings.
  double p2CA, p2CV;

};

//==========================================================================

// Helicity matrix element for the decay of a Higgs ->  two fermions.

// Because the Higgs is spin zero the Higgs production mechanism is not
// needed for calculating helicity density matrices. However, the CP mixing
// is needed.

class HMEHiggs2TwoFermions : public HelicityMatrixElement {

public:

  void initConstants();

  void initWaves(vector<HelicityParticle>&);

  complex calculateME(vector<int>);

private:

  // Coupling constants of the fermions with the Higgs.
  complex p2CA, p2CV;

};

//==========================================================================

// Base class for all tau decay helicity matrix elements.

class HMETauDecay : public HelicityMatrixElement {

public:

  virtual void initWaves(vector<HelicityParticle>&);

  virtual complex calculateME(vector<int>);

  virtual double decayWeightMax(vector<HelicityParticle>&);

protected:

  virtual void initHadronicCurrent(vector<HelicityParticle>&) {};

  virtual void calculateResonanceWeights(vector<double>&, vector<double>&,
    vector<complex>&);

};

//==========================================================================

// Helicity matrix element for a tau decaying into a single scalar meson.

class HMETau2Meson : public HMETauDecay {

public:

  void initConstants();

  void initHadronicCurrent(vector<HelicityParticle>&);

};

//==========================================================================

// Helicity matrix element for a tau decaying into two leptons.

class HMETau2TwoLeptons : public HMETauDecay {

public:

  void initConstants();

  void initWaves(vector<HelicityParticle>&);

  complex calculateME(vector<int>);

};

//==========================================================================

// Helicity matrix element for a tau decaying into two mesons through a
// vector meson resonance.

class HMETau2TwoMesonsViaVector : public HMETauDecay {

public:

  void initConstants();

  void initHadronicCurrent(vector<HelicityParticle>&);

private:

  // Resonance masses, widths, and weights.
  vector<double>  vecM, vecG, vecP, vecA;
  vector<complex> vecW;

};

//==========================================================================

// Helicity matrix element for a tau decay into two mesons through a vector
// or scalar meson resonance.

class HMETau2TwoMesonsViaVectorScalar : public HMETauDecay {

public:

  HMETau2TwoMesonsViaVectorScalar() : scaC(), vecC() {}

  void initConstants();

  void initHadronicCurrent(vector<HelicityParticle>&);

private:

  // Coupling to vector and scalar resonances.
  double scaC, vecC;

  // Resonance masses, widths, and weights.
  vector<double>  scaM, scaG, scaP, scaA, vecM, vecG, vecP, vecA;
  vector<complex> scaW, vecW;

};

//==========================================================================

// Helicity matrix element for a tau decay into three mesons (base class).

class HMETau2ThreeMesons : public HMETauDecay {

public:

  HMETau2ThreeMesons() : s1(), s2(), s3(), s4() {}

  void initConstants();

  void initHadronicCurrent(vector<HelicityParticle>&);

protected:

  // Decay mode of the tau.
  enum Mode{Pi0Pi0Pim, PimPimPip, Pi0PimK0b, PimPipKm, Pi0PimEta, PimKmKp,
            Pi0K0Km, KlPimKs, Pi0Pi0Km, KlKlPim, PimKsKs, PimK0bK0, Uknown};
  Mode mode;

  // Initialize decay mode and resonance constants (called by initConstants).
  virtual void initMode();
  virtual void initResonances() {;}

  // Initialize the momenta.
  virtual void initMomenta(vector<HelicityParticle>&);

  // Center of mass energies and momenta.
  double s1, s2, s3, s4;
  Wave4  q, q2, q3, q4;

  // Stored a1 Breit-Wigner (for speed).
  complex a1BW;

  // Form factors.
  virtual complex F1() {return complex(0, 0);}
  virtual complex F2() {return complex(0, 0);}
  virtual complex F3() {return complex(0, 0);}
  virtual complex F4() {return complex(0, 0);}

  // Phase space and Breit-Wigner for the a1.
  virtual double  a1PhaseSpace(double);
  virtual complex a1BreitWigner(double);

  // Sum running p and fixed width Breit-Wigner resonances.
  complex T(double m0, double m1, double s,
            vector<double>& M, vector<double>& G, vector<double>& W);
  complex T(double s, vector<double>& M, vector<double>& G, vector<double>& W);

};

//==========================================================================

// Helicity matrix element for a tau decay into three pions.

class HMETau2ThreePions : public HMETau2ThreeMesons {

public:

  HMETau2ThreePions() : f0M(), f0G(), f0P(), f0A(), f2M(), f2G(), f2P(),
    f2A(), sigM(), sigG(), sigP(), sigA() {}

private:

  void initResonances();

  // Resonance masses, widths, and weights.
  vector<double>  rhoM, rhoG, rhoPp, rhoAp, rhoPd, rhoAd;
  double          f0M, f0G, f0P, f0A, f2M, f2G, f2P, f2A;
  double          sigM, sigG, sigP, sigA;
  vector<complex> rhoWp, rhoWd;
  complex         f0W, f2W, sigW;

  // Form factors.
  complex F1();
  complex F2();
  complex F3();

  // Running width and Breit-Wigner for the a1.
  double  a1PhaseSpace(double);
  complex a1BreitWigner(double);

};

//==========================================================================

// Helicity matrix element for a tau decay into three mesons with kaons.

class HMETau2ThreeMesonsWithKaons : public HMETau2ThreeMesons {

public:

  HMETau2ThreeMesonsWithKaons() : kM(), piM(), piW() {}

private:

  void initResonances();

  // Resonance masses, widths, and weights.
  vector<double> rhoMa, rhoGa, rhoWa, rhoMv, rhoGv, rhoWv;
  vector<double> kstarMa, kstarGa, kstarWa, kstarMv, kstarGv, kstarWv;
  vector<double> k1Ma, k1Ga, k1Wa, k1Mb, k1Gb, k1Wb;
  vector<double> omegaM, omegaG, omegaW;
  double kM, piM, piW;

  // Form factors.
  complex F1();
  complex F2();
  complex F4();

};

//==========================================================================

// Helicity matrix element for a tau decay into generic three mesons.

class HMETau2ThreeMesonsGeneric : public HMETau2ThreeMesons {

public:

  HMETau2ThreeMesonsGeneric() : kM(), piM(), piW() {}

private:

  void initResonances();

  // Resonance masses, widths, and weights.
  vector<double> rhoMa, rhoGa, rhoWa, rhoMv, rhoGv, rhoWv;
  vector<double> kstarM, kstarG, kstarW, k1M, k1G, k1W;
  double kM, piM, piW;

  // Form factors.
  complex F1();
  complex F2();
  complex F4();

};

//==========================================================================

// Helicity matrix element for a tau decay into two pions and a photon.

class HMETau2TwoPionsGamma : public HMETauDecay {

public:

  HMETau2TwoPionsGamma() : piM() {}

  void initConstants();

  void initWaves(vector<HelicityParticle>&);

  complex calculateME(vector<int>);

protected:

  // Resonance masses, widths, and weights.
  vector<double>  rhoM, rhoG, rhoW, omegaM, omegaG, omegaW;
  double piM;

  // Form factor.
  complex F(double s, vector<double> M, vector<double> G, vector<double> W);

};

//==========================================================================

// Helicity matrix element for a tau decay into four pions.

class HMETau2FourPions : public HMETauDecay {

public:

  HMETau2FourPions() : a1M(), a1G(), rhoM(), rhoG(), sigM(), sigG(), omeM(),
    omeG(), picM(), pinM(), sigA(), sigP(), omeA(), omeP(), lambda2() {}

  void initConstants();

  void initHadronicCurrent(vector<HelicityParticle>& p);

private:

  // G-function form factors (fits).
  double G(int i, double s);

  // T-vector functions.
  Wave4 t1(Wave4&, Wave4&, Wave4&, Wave4&, Wave4&);
  Wave4 t2(Wave4&, Wave4&, Wave4&, Wave4&, Wave4&);
  Wave4 t3(Wave4&, Wave4&, Wave4&, Wave4&, Wave4&);

  // Breit-Wigner denominators for the intermediate mesons.
  complex  a1D(double s);
  complex rhoD(double s);
  complex sigD(double s);
  complex omeD(double s);

  // Form factors needed for the a1, rho, and omega.
  double  a1FormFactor(double s);
  double rhoFormFactor1(double s);
  double rhoFormFactor2(double s);
  double omeFormFactor(double s);

  // Masses and widths of the intermediate mesons.
  double a1M, a1G, rhoM, rhoG, sigM, sigG, omeM, omeG;

  // Masses for the pions (charged and neutral).
  double picM, pinM;

  // Amplitude, phases, and weights for mixing.
  double  sigA, sigP, omeA, omeP;
  complex sigW, omeW;

  // Cut-off for a1 form factor.
  double lambda2;

};

//==========================================================================

// Helicity matrix element for a tau decaying into five pions.

class HMETau2FivePions : public HMETauDecay {

public:

  HMETau2FivePions() : a1M(), a1G(), rhoM(), rhoG(), omegaM(), omegaG(),
    omegaW(), sigmaM(), sigmaG(), sigmaW() {}

  void initConstants();

  void initHadronicCurrent(vector<HelicityParticle>&);

private:

  // Hadronic currents.
  Wave4 Ja(Wave4 &q, Wave4 &q1, Wave4 &q2, Wave4 &q3, Wave4 &q4, Wave4 &q5);
  Wave4 Jb(Wave4 &q, Wave4 &q1, Wave4 &q2, Wave4 &q3, Wave4 &q4, Wave4 &q5);

  // Simplified s-wave Breit-Wigner assuming massless products.
  complex breitWigner(double s, double M, double G);

  // Masses and widths of intermediates.
  double a1M, a1G, rhoM, rhoG, omegaM, omegaG, omegaW, sigmaM, sigmaG, sigmaW;

};

//==========================================================================

// Helicity matrix element for a tau decay into flat phase space.

class HMETau2PhaseSpace : public HMETauDecay {

public:

  void initWaves(vector<HelicityParticle>&) {};

  complex calculateME(vector<int>) {return 1;}

  void calculateD(vector<HelicityParticle>&) {};

  void calculateRho(unsigned int, vector<HelicityParticle>&) {};

  double decayWeight(vector<HelicityParticle>&) {return 1.0;}

  double decayWeightMax(vector<HelicityParticle>&) {return 1.0;}

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_HelicityMatrixElements_H
