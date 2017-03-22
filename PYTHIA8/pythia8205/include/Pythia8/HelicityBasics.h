// HelicityBasics.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for a number of helper classes used in tau decays.

#ifndef Pythia8_HelicityBasics_H
#define Pythia8_HelicityBasics_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/PythiaComplex.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// The Wave4 class provides a class for complex four-vector wave functions.
// The Wave4 class can be multiplied with the GammaMatrix class to allow
// for the writing of helicity matrix elements.

class Wave4 {

public:

  // Constructors and destructor.
  Wave4() {};
  Wave4(complex v0, complex v1, complex v2, complex v3) {val[0] = v0;
    val[1] = v1; val[2] = v2; val[3] = v3;}
  Wave4(Vec4 v) {val[0] = v.e(); val[1] = v.px(); val[2] = v.py();
    val[3] = v.pz();}
  ~Wave4() {};

  // Access an element of the wave vector.
  complex& operator() (int i) {return val[i];}

  // Wave4 + Wave4.
  Wave4 operator+(Wave4 w) {return Wave4( val[0] + w.val[0],
    val[1] + w.val[1], val[2] + w.val[2], val[3] + w.val[3]);}

  // Wave4 - Wave4.
  Wave4 operator-(Wave4 w) {return Wave4( val[0] - w.val[0],
    val[1] - w.val[1], val[2] - w.val[2], val[3] - w.val[3]);}

  // - Wave4.
  Wave4 operator-() {return Wave4(-val[0], -val[1], -val[2], -val[3]);}

  // Wave4 * Wave4.
  complex operator*(Wave4 w) {return val[0] * w.val[0]
    + val[1] * w.val[1] + val[2] * w.val[2] + val[3] * w.val[3];}

  // Wave4 * complex.
  Wave4 operator*(complex s) {return Wave4(val[0] * s, val[1] * s,
    val[2] * s, val[3] * s);}

  // complex * Wave4.
  friend Wave4 operator*(complex s, const Wave4& w);

  // Wave4 * double.
  Wave4 operator*(double s) {return Wave4(val[0] * s, val[1] * s,
    val[2] * s, val[3] * s);}

  // double * Wave4.
  friend Wave4 operator*(double s, const Wave4& w);

  // Wave4 / complex.
  Wave4 operator/(complex s) {return Wave4(val[0] / s, val[1] / s,
    val[2] / s, val[3] / s);}

  // Wave4 / double.
  Wave4 operator/(double s) {return Wave4(val[0] / s, val[1] / s,
    val[2]/s, val[3]/s);}

  // Complex conjugate.
  friend Wave4 conj(Wave4 w);

  // Permutation operator.
  friend Wave4 epsilon(Wave4 w1, Wave4 w2, Wave4 w3);

  // Invariant squared mass for REAL Wave4 (to save time).
  friend double m2(Wave4 w);
  friend double m2(Wave4 w1, Wave4 w2);

  // Wave4 * GammaMatrix multiplication is defined in the GammaMatrix class.

  // Print a Wave4 vector.
  friend ostream& operator<<(ostream& output, Wave4 w);

protected:

  complex val[4];

};

//--------------------------------------------------------------------------

// Namespace function declarations; friends of Wave4 class.
Wave4 operator*(complex s, const Wave4& w);
Wave4 operator*(double s, const Wave4& w);
Wave4 conj(Wave4 w);
Wave4 epsilon(Wave4 w1, Wave4 w2, Wave4 w3);
double m2(Wave4 w);
double m2(Wave4 w1, Wave4 w2);
ostream& operator<< (ostream& os, Wave4 w);

//==========================================================================

// The GammaMatrix class is a special sparse matrix class used to write
// helicity matrix elements in conjuction with the Wave4 class. Note that
// only left to right multplication of Wave4 vectors with the GammaMatrix
// class is allowed. Additionally, subtracting a scalar from a GammaMatrix
// (or subtracting a GammaMatrix from a scalar) subtracts the scalar from
//each non-zero element of the GammaMatrix. This is designed specifically
// with the (1 - gamma^5) structure of matrix elements in mind.

class GammaMatrix {

public:

  // Constructors and destructor.
  GammaMatrix() {};
  GammaMatrix(int mu);
  ~GammaMatrix() {};

  // Access an element of the matrix.
  complex& operator() (int I, int J) {if (index[J] == I) return val[J];
    else return COMPLEXZERO; }

  // Wave4 * GammaMatrix.
  friend Wave4 operator*(Wave4 w, GammaMatrix g);

  // GammaMatrix * Scalar.
  GammaMatrix operator*(complex s) {val[0] = s*val[0]; val[1] = s*val[1];
    val[2] = s*val[2]; val[3] = s*val[3]; return *this;}

  // Scalar * GammaMatrix.
  friend GammaMatrix operator*(complex s, GammaMatrix g);

  // Gamma5 - I * Scalar.
  GammaMatrix operator-(complex s) {val[0] = val[0] - s; val[1] = val[1] - s;
    val[2] = val[2] - s; val[3] = val[3] - s; return *this;}

  // I * Scalar - Gamma5.
  friend GammaMatrix operator-(complex s, GammaMatrix g);

  // Gamma5 + I * Scalar
  GammaMatrix operator+(complex s) {val[0] = val[0] + s; val[1] = val[1] + s;
    val[2] = val[2] + s; val[3] = val[3] + s; return *this;}

  // I * Scalar + Gamma5
  friend GammaMatrix operator+(complex s, GammaMatrix g);

  // << GammaMatrix.
  friend ostream& operator<< (ostream& os, GammaMatrix g);

protected:

  complex val[4];
  int     index[4];

  // Need to define complex 0 as a variable for operator() to work.
  complex COMPLEXZERO;

};

//--------------------------------------------------------------------------

// Namespace function declarations; friends of GammaMatrix class.
Wave4 operator*(Wave4 w, GammaMatrix g);
GammaMatrix operator*(complex s, GammaMatrix g);
GammaMatrix operator-(complex s, GammaMatrix g);
GammaMatrix operator+(complex s, GammaMatrix g);
ostream& operator<< (ostream& os, GammaMatrix g);

//==========================================================================

// Helicity particle class containing helicity information, derived from
// particle base class.

class HelicityParticle : public Particle {

public:

  // Constructors.
  HelicityParticle() : Particle() { direction = 1;}
  HelicityParticle(int idIn, int statusIn = 0, int mother1In = 0,
    int mother2In = 0, int daughter1In = 0, int daughter2In = 0,
    int colIn = 0, int acolIn = 0, double pxIn = 0.,
    double pyIn = 0., double pzIn = 0., double eIn = 0.,
    double mIn = 0., double scaleIn = 0., ParticleData* ptr = 0)
    : Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In,
    colIn, acolIn, pxIn, pyIn, pzIn, eIn, mIn, scaleIn) {
    if (ptr) setPDEPtr( ptr->particleDataEntryPtr( idIn) );
    initRhoD();
    direction = 1; }
  HelicityParticle(int idIn, int statusIn, int mother1In, int mother2In,
    int daughter1In, int daughter2In, int colIn, int acolIn, Vec4 pIn,
    double mIn = 0., double scaleIn = 0., ParticleData* ptr = 0)
    : Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In,
    colIn, acolIn, pIn, mIn, scaleIn) {
    if (ptr) setPDEPtr( ptr->particleDataEntryPtr( idIn) );
    initRhoD();
    direction = 1; }
  HelicityParticle(const Particle& ptIn, ParticleData* ptr = 0)
    : Particle(ptIn) {
    indexSave = ptIn.index();
    if (ptr) setPDEPtr( ptr->particleDataEntryPtr( id()) );
    initRhoD();
    direction = 1; }

  // Methods.
  Wave4 wave(int h);
  Wave4 waveBar(int h);
  void normalize(vector< vector<complex> >& m);
  int spinStates();

  // Return and set mass (redefine from Particle).
  double m() {return mSave;}
  void   m(double mIn) {mSave = mIn; initRhoD();}

  // Event record position (redefine from Particle).
  int  index() const {return indexSave;}
  void index(int indexIn) {indexSave = indexIn;}

  // Flag for whether particle is incoming (-1) or outgoing (1).
  int direction;

  // Helicity density matrix.
  vector< vector<complex> > rho;

  // Decay matrix.
  vector< vector<complex> > D;

private:

  // Initialize the helicity density and decay matrix.
  void initRhoD() {
    rho = vector< vector<complex> >(spinStates(),
      vector<complex>(spinStates(), 0));
    D   = vector< vector<complex> >(spinStates(),
      vector<complex>(spinStates(), 0));
    for (int i = 0; i < spinStates(); i++) {
      rho[i][i] = 1.0 / spinStates(); D[i][i] = 1;}
  }

  // Particle index in the event record.
  int indexSave;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_HelicityBasics_H
