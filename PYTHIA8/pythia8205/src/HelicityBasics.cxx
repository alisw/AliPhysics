// HelicityBasics.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for helper classes
// used in tau decays.

#include "Pythia8/HelicityBasics.h"

namespace Pythia8 {

//==========================================================================

// Wave4 class.
// Friend methods to it.

//--------------------------------------------------------------------------

// complex * Wave4.

Wave4 operator*(complex s, const Wave4& w) {

  return Wave4( s * w.val[0], s * w.val[1], s * w.val[2], s * w.val[3]);

}

//--------------------------------------------------------------------------

// double * Wave4.

Wave4 operator*(double s, const Wave4& w) {

  return Wave4( s * w.val[0], s * w.val[1], s * w.val[2], s * w.val[3]);

}

//--------------------------------------------------------------------------

// Complex conjugate.

Wave4 conj(Wave4 w) {

  w(0) = conj(w(0));
  w(1) = conj(w(1));
  w(2) = conj(w(2));
  w(3) = conj(w(3));
  return w;

}

//--------------------------------------------------------------------------

// Permutation operator.

Wave4 epsilon(Wave4 w1, Wave4 w2, Wave4 w3) {

  Wave4 w4;
  w4(0) = -(w1(1) * w2(2) * w3(3)) + (w1(1) * w2(3) * w3(2))
    + (w1(2) * w2(1) * w3(3)) - (w1(2) * w2(3) * w3(1))
    - (w1(3) * w2(1) * w3(2)) + (w1(3) * w2(2) * w3(1));
  w4(1) = -(w1(0) * w2(2) * w3(3)) + (w1(0) * w2(3) * w3(2))
    + (w1(2) * w2(0) * w3(3)) - (w1(2) * w2(3) * w3(0))
    - (w1(3) * w2(0) * w3(2)) + (w1(3) * w2(2) * w3(0));
  w4(2) = (w1(0) * w2(1) * w3(3)) - (w1(0) * w2(3) * w3(1))
    - (w1(1) * w2(0) * w3(3)) + (w1(1) * w2(3) * w3(0))
    + (w1(3) * w2(0) * w3(1)) - (w1(3) * w2(1) * w3(0));
  w4(3) = -(w1(0) * w2(1) * w3(2)) + (w1(0) * w2(2) * w3(1))
    + (w1(1) * w2(0) * w3(2)) - (w1(1) * w2(2) * w3(0))
    - (w1(2) * w2(0) * w3(1)) + (w1(2) * w2(1) * w3(0));
  return w4;

}

//--------------------------------------------------------------------------

// Invariant squared mass for REAL Wave4 (to save time).

double m2(Wave4 w) {

  return real(w(0)) * real(w(0)) - real(w(1)) * real(w(1))
    - real(w(2)) * real(w(2)) - real(w(3)) * real(w(3));

}

double m2(Wave4 w1, Wave4 w2) {

  return real(w1(0)) * real(w2(0)) - real(w1(1)) * real(w2(1))
       - real(w1(2)) * real(w2(2)) - real(w1(3)) * real(w2(3));

}

//--------------------------------------------------------------------------

// Print a Wave4 vector.

ostream& operator<< (ostream& os, Wave4 w) {

  os << left << setprecision(2);
  for (int i = 0; i < 4; i++) os << setw(20) << w.val[i];
  os << "\n";
  return os;

}

//==========================================================================

// Constructor for the GammaMatrix class. Gamma(1) through Gamma(3) give the
// corresponding gamma matrices using the Weyl basis as outlined in the HELAS
// paper. Gamma(4) gives the +--- metric, while Gamma(5) gives the gamma^5
// matrix.

GammaMatrix::GammaMatrix(int mu) {

  COMPLEXZERO = complex( 0., 0.);

  if (mu == 0) {
    val[0] =  1.; val[1] =  1.; val[2] =  1.; val[3] =  1.;
    index[0] = 2; index[1] = 3; index[2] = 0; index[3] = 1;

  } else if (mu == 1) {
    val[0] = -1.; val[1] = -1.; val[2]  = 1.; val[3]  = 1.;
    index[0] = 3; index[1] = 2; index[2] = 1; index[3] = 0;

  } else if (mu == 2) {
    val[0] = complex(0.,-1.); val[1] = complex(0.,1.);
    val[2] = complex(0.,1.);  val[3] = complex(0.,-1.);
    index[0] = 3; index[1] = 2; index[2] = 1; index[3] = 0;

  } else if (mu == 3) {
    val[0] = -1.; val[1] =  1.; val[2] =  1.; val[3] = -1.;
    index[0] = 2; index[1] = 3; index[2] = 0; index[3] = 1;

  } else if (mu == 4) {
    val[0] =  1.; val[1] = -1.; val[2] = -1.; val[3] = -1.;
    index[0] = 0; index[1] = 1; index[2] = 2; index[3] = 3;

  } else if (mu == 5) {
    val[0] = -1.; val[1] = -1.; val[2] =  1.; val[3] =  1.;
    index[0] = 0; index[1] = 1; index[2] = 2; index[3] = 3;
  }

}

//--------------------------------------------------------------------------

// Wave4 * GammaMatrix.

Wave4 operator*(Wave4 w, GammaMatrix g) {

  complex w0 = w(g.index[0]);
  complex w1 = w(g.index[1]);
  complex w2 = w(g.index[2]);
  complex w3 = w(g.index[3]);
  w(0) = w0 * g.val[0];
  w(1) = w1 * g.val[1];
  w(2) = w2 * g.val[2];
  w(3) = w3 * g.val[3];
  return w;

}

//--------------------------------------------------------------------------

// Scalar * GammaMatrix.

GammaMatrix operator*(complex s, GammaMatrix g) {

  g.val[0] = s * g.val[0];
  g.val[1] = s*g.val[1];
  g.val[2] = s * g.val[2];
  g.val[3] = s*g.val[3];
  return g;

}

//--------------------------------------------------------------------------

// I * Scalar - Gamma5.

GammaMatrix operator-(complex s, GammaMatrix g) {

  g.val[0] = s - g.val[0];
  g.val[1] = s - g.val[1];
  g.val[2] = s - g.val[2];
  g.val[3] = s - g.val[3];
  return g;

}

//--------------------------------------------------------------------------

// I * Scalar + Gamma5.

GammaMatrix operator+(complex s, GammaMatrix g) {

  g.val[0] = s + g.val[0];
  g.val[1] = s + g.val[1];
  g.val[2] = s + g.val[2];
  g.val[3] = s + g.val[3];
  return g;

}

//--------------------------------------------------------------------------

// Print GammaMatrix.

ostream& operator<< (ostream& os, GammaMatrix g) {

  os << left << setprecision(2);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) os << setw(20) << g(i,j);
    os << "\n";
  }
  return os;

}

//==========================================================================

// Weyl helicity wave functions for spin 1/2 fermions and spin 1 boson.

// This is the basis as given by the HELAS collaboration on page 122 of
// "HELas: HELicity Amplitude Subroutines for Feynman Diagram Evaluations"
// by H. Murayama, I. Watanabe, K. Hagiwara.

// The spinors become ill-defined for p -> -pz and the polarization vectors
// become ill-defined when pT -> 0. For these special cases limits are used.

//--------------------------------------------------------------------------

// Return wave vector for given helicity.

Wave4 HelicityParticle::wave(int h) {

  // Create wave vector to return.
  Wave4 w;

  // Fermion (spin 1/2) spinor.
  if (spinType() == 2) {

    // Calculate helicity independent normalization.
    double P     = pAbs();
    double n     = sqrtpos(2*P*(P+pz()));
    bool aligned = P + pz() == 0;

    // Calculate eigenspinor basis.
    vector< vector<complex> > xi(2, vector<complex>(2));
    // Helicity -1
    xi[0][0] = aligned ? -1 : complex(-px(),py())/n;
    xi[0][1] = aligned ?  0 : (P+pz())/n;
    // Helicity +1
    xi[1][0] = aligned ? 0 : (P+pz())/n;
    xi[1][1] = aligned ? 1 : complex(px(),py())/n;

    // Calculate helicity dependent normalization.
    vector<double> omega(2);
    omega[0] = sqrtpos(e()-P);
    omega[1] = sqrtpos(e()+P);
    vector<double> hsign(2,1);
    hsign[0] = -1;

    // Create particle spinor.
    if (this->id() > 0) {
      w(0) = omega[!h] * xi[h][0];
      w(1) = omega[!h] * xi[h][1];
      w(2) = omega[h]  * xi[h][0];
      w(3) = omega[h]  * xi[h][1];

    // Create anti-particle spinor.
    } else {
      w(0) = hsign[!h] * omega[h]  * xi[!h][0];
      w(1) = hsign[!h] * omega[h]  * xi[!h][1];
      w(2) = hsign[h]  * omega[!h] * xi[!h][0];
      w(3) = hsign[h]  * omega[!h] * xi[!h][1];
    }

  // Boson (spin 1) polarization vector.
  } else if (spinType() == 3) {
    double P  = pAbs();
    double PT = pT();

    // Create helicity +1 or -1 polarization vector.
    if (h >= 0 && h <= 1) {
      double hsign = h ? -1 : 1;
      if (P == 0) {
        w(0) = 0;
        w(1) = hsign / sqrt(2);
        w(2) = complex(0, 1/sqrt(2));
        w(3) = 0;
      } else if (PT == 0) {
        w(0) = 0;
        w(1) = hsign / sqrt(2);
        w(2) = complex(0, (pz() > 0 ? 1 : -1) / sqrt(2));
        w(3) = complex(-hsign * PT / P, 0) / sqrt(2);
      } else {
        w(0) = 0;
        w(1) = complex(hsign * px() * pz() / (P * PT), -py() / PT) / sqrt(2);
        w(2) = complex(hsign * py() * pz() / (P * PT),  px() / PT) / sqrt(2);
        w(3) = complex(-hsign * PT / P, 0) / sqrt(2);
      }

    // Create helicity 0 polarization vector (ensure boson massive).
    } else if (h == 2 && spinStates() == 3) {
      if (P == 0) {
        w(0) = 0;
        w(1) = 0;
        w(2) = 0;
        w(3) = 1;
      } else {
        w(0) = P / m();
        w(1) = px() * e() / (m() * P);
        w(2) = py() * e() / (m() * P);
        w(3) = pz() * e() / (m() * P);
      }
    }

  // Unknown wave function.
  } else {
    w(0) = 0;
    w(1) = 0;
    w(2) = 0;
    w(3) = 0;
  }

  // Done.
  return w;

}

//--------------------------------------------------------------------------

// Bar of a wave function.

Wave4 HelicityParticle::waveBar(int h) {

  if (spinType() == 2) return conj(wave(h)) * GammaMatrix(0);
  else                 return conj(wave(h));

}

//--------------------------------------------------------------------------

// Normalize the rho or D matrices.

void HelicityParticle::normalize(vector< vector<complex> >& matrix) {

  complex trace = 0;
  for (unsigned int i = 0; i < matrix.size(); i++) trace += matrix[i][i];
  for (unsigned int i = 0; i < matrix.size(); i++) {
    for (unsigned int j = 0; j < matrix.size(); j++) {
      if (trace != complex(0,0)) matrix[i][j] /= trace;
      else matrix[i][j] = 1 / static_cast<double>(matrix.size());
    }
  }

}

//--------------------------------------------------------------------------

// Return the number of spin states.

  int HelicityParticle::spinStates() {

    int sT = spinType();
    if (sT == 0) return 1;
    else if (sT != 2 && m() == 0) return sT - 1;
    else return sT;

  }

//==========================================================================

} // end namespace Pythia8
