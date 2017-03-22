// HelicityMatrixElements.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for physics classes
// used in tau decays.

#include "Pythia8/HelicityMatrixElements.h"

namespace Pythia8 {

//==========================================================================

// The HelicityMatrixElements class.

//--------------------------------------------------------------------------

// Initialize the helicity matrix element.

void HelicityMatrixElement::initPointers(ParticleData* particleDataPtrIn,
  Couplings* couplingsPtrIn, Settings* settingsPtrIn) {

  particleDataPtr = particleDataPtrIn;
  couplingsPtr    = couplingsPtrIn;
  settingsPtr     = settingsPtrIn;
  for(int i = 0; i <= 5; i++)
    gamma.push_back(GammaMatrix(i));

}

//--------------------------------------------------------------------------

// Initialize the channel for the helicity matrix element.

HelicityMatrixElement* HelicityMatrixElement::initChannel(
  vector<HelicityParticle>& p) {

  pID.clear();
  pM.clear();
  for(int i = 0; i < static_cast<int>(p.size()); i++) {
    pID.push_back(p[i].id());
    pM.push_back(p[i].m());
  }
  initConstants();
  return this;

}

//--------------------------------------------------------------------------

// Calculate a particle's decay matrix.

void HelicityMatrixElement::calculateD(vector<HelicityParticle>& p) {

  // Reset the D matrix to zero.
  for (int i = 0; i < p[0].spinStates(); i++) {
    for (int j = 0; j < p[0].spinStates(); j++) {
        p[0].D[i][j] = 0;
    }
  }

  // Initialize the wave functions.
  initWaves(p);

  // Create the helicity vectors.
  vector<int> h1(p.size(),0);
  vector<int> h2(p.size(),0);

  // Call the recursive sub-method.
  calculateD(p, h1, h2, 0);

  // Normalize the decay matrix.
  p[0].normalize(p[0].D);

}

//--------------------------------------------------------------------------

// Recursive sub-method for calculating a particle's decay matrix.

void HelicityMatrixElement::calculateD(vector<HelicityParticle>& p,
  vector<int>& h1, vector<int>& h2, unsigned int i) {

  if (i < p.size()) {
    for (h1[i] = 0; h1[i] < p[i].spinStates(); h1[i]++) {
        for (h2[i] = 0; h2[i] < p[i].spinStates(); h2[i]++) {
          calculateD(p, h1, h2, i+1);
        }
    }
  }
  else {
    p[0].D[h1[0]][h2[0]] += calculateME(h1) * conj(calculateME(h2)) *
        calculateProductD(p, h1, h2);
  }

}

//--------------------------------------------------------------------------

// Calculate a particle's helicity density matrix.

void HelicityMatrixElement::calculateRho(unsigned int idx,
  vector<HelicityParticle>& p) {

  // Reset the rho matrix to zero.
  for (int i = 0; i < p[idx].spinStates(); i++) {
    for (int j = 0; j < p[idx].spinStates(); j++) {
        p[idx].rho[i][j] = 0;
    }
  }

  // Initialize the wave functions.
  initWaves(p);

  // Create the helicity vectors.
  vector<int> h1(p.size(),0);
  vector<int> h2(p.size(),0);

  // Call the recursive sub-method.
  calculateRho(idx, p, h1, h2, 0);

  // Normalize the density matrix.
  p[idx].normalize(p[idx].rho);

}

//--------------------------------------------------------------------------

// Recursive sub-method for calculating a particle's helicity density matrix.

void HelicityMatrixElement::calculateRho(unsigned int idx,
  vector<HelicityParticle>& p, vector<int>& h1, vector<int>& h2,
  unsigned int i) {

  if (i < p.size()) {
    for (h1[i] = 0; h1[i] < p[i].spinStates(); h1[i]++) {
        for (h2[i] = 0; h2[i] < p[i].spinStates(); h2[i]++) {
          calculateRho(idx, p, h1, h2, i+1);
        }
    }
  }
  else {
    // Calculate rho from a hard process.
    if (p[1].direction < 0)
        p[idx].rho[h1[idx]][h2[idx]] += p[0].rho[h1[0]][h2[0]] *
          p[1].rho[h1[1]][h2[1]] * calculateME(h1)*conj(calculateME(h2)) *
          calculateProductD(idx, 2, p, h1, h2);
    // Calculate rho from a decay.
    else
        p[idx].rho[h1[idx]][h2[idx]] += p[0].rho[h1[0]][h2[0]] *
          calculateME(h1)*conj(calculateME(h2)) *
          calculateProductD(idx, 1, p, h1, h2);
    return;
  }

}

//--------------------------------------------------------------------------

// Calculate a decay's weight.

double HelicityMatrixElement::decayWeight(vector<HelicityParticle>& p) {

  complex weight = complex(0,0);

  // Initialize the wave functions.
  initWaves(p);

  // Create the helicity vectors.
  vector<int> h1(p.size(),0);
  vector<int> h2(p.size(),0);

  // Call the recursive sub-method.
  decayWeight(p, h1, h2, weight, 0);

  return real(weight);

}

//--------------------------------------------------------------------------

// Recursive sub-method for calculating a decay's weight.

void HelicityMatrixElement::decayWeight(vector<HelicityParticle>& p,
  vector<int>& h1, vector<int>& h2, complex& weight, unsigned int i) {

  if (i < p.size()) {
    for (h1[i] = 0; h1[i] < p[i].spinStates(); h1[i]++) {
        for (h2[i] = 0; h2[i] < p[i].spinStates(); h2[i]++) {
          decayWeight(p, h1, h2, weight, i+1);
        }
    }
  }
  else {
    weight += p[0].rho[h1[0]][h2[0]] * calculateME(h1) *
        conj(calculateME(h2)) * calculateProductD(p, h1, h2);
  }

}

//--------------------------------------------------------------------------

// Calculate the product of the decay matrices (hard process).

complex HelicityMatrixElement::calculateProductD(unsigned int idx,
  unsigned int start, vector<HelicityParticle>& p,
  vector<int>& h1, vector<int>& h2) {

  complex answer(1,0);
  for (unsigned int i = start; i < p.size(); i++) {
    if (i != idx) {
        answer *= p[i].D[h1[i]][h2[i]];
    }
  }
  return answer;

}

//--------------------------------------------------------------------------

// Calculate the product of the decay matrices (decay process).

complex HelicityMatrixElement::calculateProductD(
  vector<HelicityParticle>& p, vector<int>& h1, vector<int>& h2) {

  complex answer(1,0);
  for (unsigned int i = 1; i < p.size(); i++) {
    answer *= p[i].D[h1[i]][h2[i]];
  }
  return answer;

}

//--------------------------------------------------------------------------

// Initialize a fermion line.

void HelicityMatrixElement::setFermionLine(int position,
  HelicityParticle& p0, HelicityParticle& p1) {

  vector< Wave4 > u0, u1;

  // First particle is incoming and particle, or outgoing and anti-particle.
  if (p0.id()*p0.direction < 0) {
    pMap[position] = position; pMap[position+1] = position+1;
    for (int h = 0; h < p0.spinStates(); h++) u0.push_back(p0.wave(h));
    for (int h = 0; h < p1.spinStates(); h++) u1.push_back(p1.waveBar(h));
  }
  // First particle is outgoing and particle, or incoming and anti-particle.
  else {
    pMap[position] = position+1; pMap[position+1] = position;
    for (int h = 0; h < p0.spinStates(); h++) u1.push_back(p0.waveBar(h));
    for (int h = 0; h < p1.spinStates(); h++) u0.push_back(p1.wave(h));
  }
  u.push_back(u0); u.push_back(u1);

}

//--------------------------------------------------------------------------

// Return a fixed width Breit-Wigner.

complex HelicityMatrixElement::breitWigner(double s, double M, double G) {

  return (-M*M + complex(0, 1) * M * G) / (s - M*M + complex(0, 1) * M * G);

}

//--------------------------------------------------------------------------

// Return an s-wave BreitWigner.

complex HelicityMatrixElement::sBreitWigner(double m0, double m1, double s,
   double M, double G) {

  double gs = sqrtpos((s - pow2(m0+m1)) * (s - pow2(m0-m1))) / (2*sqrtpos(s));
  double gM = sqrtpos((M*M - pow2(m0+m1)) * (M*M - pow2(m0-m1))) / (2*M);
  return M*M / (M*M - s - complex(0,1)*G*M*M/sqrtpos(s)*(gs/gM));

}

//--------------------------------------------------------------------------

// Return a p-wave BreitWigner.

complex HelicityMatrixElement::pBreitWigner(double m0, double m1, double s,
  double M, double G) {

  double gs = sqrtpos((s - pow2(m0+m1)) * (s - pow2(m0-m1))) / (2*sqrtpos(s));
  double gM = sqrtpos((M*M - pow2(m0+m1)) * (M*M - pow2(m0-m1))) / (2*M);
  return M*M / (M*M - s - complex(0,1)*G*M*M/sqrtpos(s)*pow3(gs/gM));

}

//--------------------------------------------------------------------------

// Return a d-wave BreitWigner.

complex HelicityMatrixElement::dBreitWigner(double m0, double m1, double s,
  double M, double G) {

  double gs = sqrtpos((s - pow2(m0+m1)) * (s - pow2(m0-m1))) / (2*sqrtpos(s));
  double gM = sqrtpos((M*M - pow2(m0+m1)) * (M*M - pow2(m0-m1))) / (2*M);
  return M*M / (M*M - s - complex(0,1)*G*M*M/sqrtpos(s)*pow5(gs/gM));

}

//==========================================================================

// Helicity matrix element for two fermions -> W -> two fermions. This matrix
// element handles s-channel hard processes in addition to t-channel, assuming
// the first two particles are a fermion line and the second two particles
// are a fermion line. This matrix element is not scaled with respect to W
// propagator energy as currently this matrix element is used only for
// calculating helicity density matrices.

//--------------------------------------------------------------------------

// Initialize spinors for the helicity matrix element.

void HMETwoFermions2W2TwoFermions::initWaves(vector<HelicityParticle>& p) {

  u.clear();
  pMap.resize(4);
  setFermionLine(0,p[0],p[1]);
  setFermionLine(2,p[2],p[3]);

}

//--------------------------------------------------------------------------

  // Return element for the helicity matrix element.

complex HMETwoFermions2W2TwoFermions::calculateME(vector<int> h) {

  complex answer(0,0);
  for (int mu = 0; mu <= 3; mu++) {
    answer += (u[1][h[pMap[1]]] * gamma[mu] * (1 - gamma[5])
      * u[0][h[pMap[0]]]) * gamma[4](mu,mu) * (u[3][h[pMap[3]]]
      * gamma[mu] * (1 - gamma[5]) * u[2][h[pMap[2]]]);
  }
  return answer;

}

//==========================================================================

// Helicity matrix element for two fermions -> photon -> two fermions. This
// matrix element can be combined with the Z matrix element to provide full
// interference effects.

// p0Q: charge of the incoming fermion line
// p2Q: charge of the outgoing fermion line
// s: center of mass energy

//--------------------------------------------------------------------------

// Initialize wave functions for the helicity matrix element.

void HMETwoFermions2Gamma2TwoFermions::initWaves(
  vector<HelicityParticle>& p) {

  u.clear();
  pMap.resize(4);
  setFermionLine(0, p[0], p[1]);
  setFermionLine(2, p[2], p[3]);
  s = max( 1., pow2(p[4].m()));
  p0Q = p[0].charge(); p2Q = p[2].charge();

}

//--------------------------------------------------------------------------

// Return element for the helicity matrix element.


complex HMETwoFermions2Gamma2TwoFermions::calculateME(vector<int> h) {

  complex answer(0,0);
  for (int mu = 0; mu <= 3; mu++) {
    answer += (u[1][h[pMap[1]]] * gamma[mu] * u[0][h[pMap[0]]])
      * gamma[4](mu,mu) * (u[3][h[pMap[3]]] * gamma[mu] * u[2][h[pMap[2]]]);
  }
  return p0Q*p2Q * answer / s;

}

//==========================================================================

// Helicity matrix element for two fermions -> Z -> two fermions. This matrix
// element can be combined with the photon matrix element to provide full
// interference effects.

// Note that there is a double contraction in the Z matrix element, which can
// be very time consuming. If the two incoming fermions are oriented along
// the z-axis, their helicities must be opposite for a non-zero matrix element
// term. Consequently, this check is made to help speed up the matrix element.

// sin2W: sine of the Weinberg angle
// cos2W: cosine of the Weinberg angle
// zM: on-shell mass of the Z
// zG: on-shell width of the Z
// p0CA: axial coupling of particle 0 to the Z
// p2CA: axial coupling of particle 2 to the Z
// p0CV: vector coupling of particle 0 to the Z
// p2CV: vector coupling of particle 2 to the Z
// zaxis: true if the incoming fermions are oriented along the z-axis

//--------------------------------------------------------------------------

// Initialize the constant for the helicity matrix element.

void HMETwoFermions2Z2TwoFermions::initConstants() {

  // Set the Weinberg angle.
  sin2W = couplingsPtr->sin2thetaW();
  cos2W = couplingsPtr->cos2thetaW();
  // Set the on-shell Z mass and width.
  zG = particleDataPtr->mWidth(23);
  zM = particleDataPtr->m0(23);
  // Set the vector and axial couplings to the fermions.
  p0CA = couplingsPtr->af(abs(pID[0]));
  p2CA = couplingsPtr->af(abs(pID[2]));
  p0CV = couplingsPtr->vf(abs(pID[0]));
  p2CV = couplingsPtr->vf(abs(pID[2]));

}

//--------------------------------------------------------------------------

// Initialize wave functions for the helicity matrix element.

void HMETwoFermions2Z2TwoFermions::initWaves(vector<HelicityParticle>& p) {

  vector< Wave4 > u4;
  u.clear();
  pMap.resize(4);
  setFermionLine(0, p[0], p[1]);
  setFermionLine(2, p[2], p[3]);
  u4.push_back(Wave4(p[2].p() + p[3].p()));
  u.push_back(u4);
  // Center of mass energy.
  s = max( 1., pow2(p[4].m()));
  // Check if incoming fermions are oriented along z-axis.
  zaxis = (p[0].pAbs() == fabs(p[0].pz())) &&
    (p[1].pAbs() == fabs(p[1].pz()));

}

//--------------------------------------------------------------------------

// Return element for helicity matrix element.

complex HMETwoFermions2Z2TwoFermions::calculateME(vector<int> h) {

  complex answer(0,0);
  // Return zero if correct helicity conditions.
  if (h[0] == h[1] && zaxis) return answer;
  for (int mu = 0; mu <= 3; mu++) {
    for (int nu = 0; nu <= 3; nu++) {
        answer +=
          (u[1][h[pMap[1]]] * gamma[mu] * (p0CV - p0CA * gamma[5]) *
           u[0][h[pMap[0]]]) *
          (gamma[4](mu,nu) - gamma[4](mu,mu)*u[4][0](mu) *
           gamma[4](nu,nu) * u[4][0](nu) / (zM*zM)) *
          (u[3][h[pMap[3]]] * gamma[nu] * (p2CV - p2CA * gamma[5]) *
           u[2][h[pMap[2]]]);
    }
  }
  return answer / (16 * pow2(sin2W * cos2W) *
                     (s - zM*zM + complex(0, s*zG/zM)));

}

//==========================================================================

// Helicity matrix element for two fermions -> photon/Z -> two fermions. Full
// interference is obtained by combining the photon and Z helicity matrix
// elements.

// In general the initPointers and initChannel methods should not be
// redeclared.

//--------------------------------------------------------------------------

// Initialize the matrix element.

void HMETwoFermions2GammaZ2TwoFermions::initPointers(
  ParticleData* particleDataPtrIn, Couplings* couplingsPtrIn,
  Settings *settingsPtrIn) {

  zHME.initPointers(particleDataPtrIn, couplingsPtrIn, settingsPtrIn);
  gHME.initPointers(particleDataPtrIn, couplingsPtrIn, settingsPtrIn);

}

//--------------------------------------------------------------------------

// Initialize the channel for the helicity matrix element.

  HelicityMatrixElement* HMETwoFermions2GammaZ2TwoFermions::initChannel(
    vector<HelicityParticle>& p) {

    zHME.initChannel(p);
    zHME.initChannel(p);
    return this;

}

//--------------------------------------------------------------------------

// Initialize wave functions for the helicity matrix element.

void HMETwoFermions2GammaZ2TwoFermions::initWaves(
  vector<HelicityParticle>& p) {

  zHME.initWaves(p);
  gHME.initWaves(p);

}

//--------------------------------------------------------------------------

// Return element for the helicity matrix element.

complex HMETwoFermions2GammaZ2TwoFermions::calculateME(vector<int> h) {

  return zHME.calculateME(h) + gHME.calculateME(h);

}

//==========================================================================

// Helicity matrix element for X -> two fermions.

// Base class for the W, photon, and Z -> two fermions helicity matrix
// elements.

//--------------------------------------------------------------------------

// Initialize wave functions for the helicity matrix element.

void HMEX2TwoFermions::initWaves(vector<HelicityParticle>& p) {

  u.clear();
  pMap.resize(4);
  // Initialize W wave function.
  vector< Wave4 > u1;
  pMap[1] = 1;
  for (int h = 0; h < p[pMap[1]].spinStates(); h++)
    u1.push_back(p[pMap[1]].wave(h));
  u.push_back(u1);
  // Initialize fermion wave functions.
  setFermionLine(2, p[2], p[3]);

}

//==========================================================================

// Helicity matrix element for W -> two fermions.

// This matrix element is used when the production of the W is from an
// unknown process.

//--------------------------------------------------------------------------

// Return element for helicity matrix element.

complex HMEW2TwoFermions::calculateME(vector<int> h) {

  complex answer(0,0);
  for (int mu = 0; mu <= 3; mu++) {
    answer +=
      u[0][h[pMap[1]]](mu) * (u[2][h[pMap[3]]] * gamma[mu]
                              * (1 - gamma[5]) *  u[1][h[pMap[2]]]);
  }
  return answer;
}

//==========================================================================

// Helicity matrix element for photon -> two fermions.

// This matrix element is used when the production of the photon is from an
// unknown process.

//--------------------------------------------------------------------------

// Return element for helicity matrix element.

complex HMEGamma2TwoFermions::calculateME(vector<int> h) {

  complex answer(0,0);
  for (int mu = 0; mu <= 3; mu++) {
    answer +=
      u[0][h[pMap[1]]](mu) * (u[2][h[pMap[3]]] * gamma[mu] * u[1][h[pMap[2]]]);
  }
  return answer;
}

//==========================================================================

// Helicity matrix element for Z -> two fermions.

// This matrix element is used when the production of the Z is from an
// unknown process.

// p2CA: axial coupling of particle 2 to the Z
// p2CV: vector coupling of particle 2 to the Z

//--------------------------------------------------------------------------

// Initialize the constant for the helicity matrix element.

void HMEZ2TwoFermions::initConstants() {

  // Set the vector and axial couplings to the fermions.
  p2CA = couplingsPtr->af(abs(pID[2]));
  p2CV = couplingsPtr->vf(abs(pID[2]));

}

//--------------------------------------------------------------------------

// Return element for helicity matrix element.

complex HMEZ2TwoFermions::calculateME(vector<int> h) {

  complex answer(0,0);
  for (int mu = 0; mu <= 3; mu++) {
    answer +=
      u[0][h[pMap[1]]](mu) * (u[2][h[pMap[3]]] * gamma[mu]
                              * (p2CV - p2CA * gamma[5]) *  u[1][h[pMap[2]]]);
  }
  return answer;
}

//==========================================================================

// Helicity matrix element for the decay of a Higgs to two fermions.

// All SM and MSSM Higgses couple to fermions with a vertex factor of
// (pfCV - pfCA * gamma[5]) where pf indicates the type of fermion line. For
// simplicity for the SM and MSSM CP even Higgses pfCV is set to one, and
// pfCA to zero, as this matrix element is used only for calculating helicity
// density matrices.

// p2CA: in the SM/MSSM this coupling is zero for CP-even and for CP-odd:
//     -g_w * m_f / (2 * m_W)
//                            * cot(beta) for A^0 u-type
//                            * tan(beta) for A^0 d-type
// p2CA: for the charged Higgs in the MSSM this coupling is given by:
//     +/- i * g / (sqrt(8.) * m_W) * (m_d * tan(beta) - m_u * cot(beta))
// p2CV: in the SM/MSSM this coupling is zero for CP-odd and for CP-even:
//     i * g_w * m_f / (2 * m_W)
//                               * -1 for the SM H
//                               * -sin(alpha) / sin(beta) for H^0 u-type
//                               * -cos(alpha) / cos(beta) for H^0 d-type
//                               * -cos(alpha) / sin(beta) for h^0 u-type
//                               *  sin(alpha) / cos(beta) for h^0 d-type
// p2CV: for the charged Higgs in the MSSM this coupling is given by:
//     i * g / (sqrt(8.) * m_W) * (m_d * tan(beta) + m_u * cot(beta))

//--------------------------------------------------------------------------

// Initialize the constant for the helicity matrix element.

void HMEHiggs2TwoFermions::initConstants() {

  // Neutral constants; settings available.
  if (settingsPtr) {
    // Set the H1 mixing.
    int mode   = settingsPtr->mode("HiggsH1:parity");
    double eta = settingsPtr->parm("HiggsH1:etaParity");
    double phi = settingsPtr->parm("HiggsH1:phiParity");
    if      (mode == 2) {p2CAH1 = 1;        p2CVH1 = 0;}
    else if (mode == 3) {p2CAH1 = eta;      p2CVH1 = 1;}
    else if (mode == 4) {p2CAH1 = cos(phi); p2CVH1 = sin(phi);}
    else                {p2CAH1 = 0;        p2CVH1 = 1;}
    // Set the H2 mixing.
    mode = settingsPtr->mode("HiggsH2:parity");
    eta  = settingsPtr->parm("HiggsH2:etaParity");
    phi  = settingsPtr->parm("HiggsH2:phiParity");
    if      (mode == 2) {p2CAH2 = 1;        p2CVH2 = 0;}
    else if (mode == 3) {p2CAH2 = eta;      p2CVH2 = 1;}
    else if (mode == 4) {p2CAH2 = cos(phi); p2CVH2 = sin(phi);}
    else                {p2CAH2 = 0;        p2CVH2 = 1;}
    // Set the A3 mixing.
    mode = settingsPtr->mode("HiggsA3:parity");
    eta  = settingsPtr->parm("HiggsA3:etaParity");
    phi  = settingsPtr->parm("HiggsA3:phiParity");
    if      (mode == 1) {p2CAA3 = 0;        p2CVA3 = 1;}
    else if (mode == 3) {p2CAA3 = eta;      p2CVA3 = 1;}
    else if (mode == 4) {p2CAA3 = cos(phi); p2CVA3 = sin(phi);}
    else                {p2CAA3 = 1;        p2CVA3 = 0;}

  // Neutral constants; default SM/MSSM.
  } else {
    p2CAH1 = 0; p2CVH1 = 1;
    p2CAH2 = 0; p2CVH2 = 1;
    p2CAA3 = 1; p2CVA3 = 0;
  }

  // Set the H4 mixing.
  p2CAH4 = 1; p2CVH4 = 1;

}

//--------------------------------------------------------------------------

// Initialize wave functions for the helicity matrix element.

void HMEHiggs2TwoFermions::initWaves(vector<HelicityParticle>& p) {

  u.clear();
  pMap.resize(4);
  if      (abs(p[1].id()) == 25) {p2CA = p2CAH1; 
                                  p2CV = complex(0, 1) * p2CVH1;}
  else if (abs(p[1].id()) == 35) {p2CA = p2CAH2; 
                                  p2CV = complex(0, 1) * p2CVH2;}
  else if (abs(p[1].id()) == 36) {p2CA = p2CAA3; 
                                  p2CV = complex(0, 1) * p2CVA3;}
  else if (abs(p[1].id()) == 37) {
    p2CA = p[1].id() == 37 ? p2CAH4 : -p2CAH4; p2CV = p2CVH4;
  } else {p2CA = 0; p2CV = 0;}
  setFermionLine(2, p[2], p[3]);

}

//--------------------------------------------------------------------------

// Return element for the helicity matrix element.

complex HMEHiggs2TwoFermions::calculateME(vector<int> h) {

  return (u[1][h[pMap[3]]] * (p2CV + p2CA * gamma[5]) * u[0][h[pMap[2]]]);

}

//==========================================================================

// Base class for all tau decay matrix elements. This class derives from
// the HelicityMatrixElement class and redefines some of the virtual functions.

// One new method, initHadronicCurrent is defined which initializes the
// hadronic current in the initWaves method. For each tau decay matrix element
// the hadronic current method must be redefined accordingly, but initWaves
// should not be redefined.

//--------------------------------------------------------------------------

// Initialize wave functions for the helicity matrix element.
void HMETauDecay::initWaves(vector<HelicityParticle>& p) {

  u.clear();
  pMap.resize(p.size());
  setFermionLine(0, p[0], p[1]);
  initHadronicCurrent(p);

}

//--------------------------------------------------------------------------

// Return element for the helicity matrix element.
complex HMETauDecay::calculateME(vector<int> h) {

  complex answer(0,0);
  for (int mu = 0; mu <= 3; mu++) {
    answer +=
        (u[1][h[pMap[1]]] * gamma[mu] * (1 - gamma[5]) * u[0][h[pMap[0]]])
        * gamma[4](mu,mu) * u[2][0](mu);
  }
  return answer;

}

//--------------------------------------------------------------------------

// Return the maximum decay weight for the helicity matrix element.

double HMETauDecay::decayWeightMax(vector<HelicityParticle>& p) {

  // Determine the maximum on-diagonal element of rho.
  double on  = real(p[0].rho[0][0]) > real(p[0].rho[1][1]) ?
    real(p[0].rho[0][0]) : real(p[0].rho[1][1]);
  // Determine the maximum off-diagonal element of rho.
  double off = fabs(real(p[0].rho[0][1])) + fabs(imag(p[0].rho[0][1]));
  return  DECAYWEIGHTMAX * (on + off);

}

//--------------------------------------------------------------------------

// Calculate complex resonance weights given a phase and amplitude vector.

void HMETauDecay::calculateResonanceWeights(vector<double>& phase,
  vector<double>& amplitude, vector<complex>& weight) {

  for (unsigned int i = 0; i < phase.size(); i++)
    weight.push_back(amplitude[i] * (cos(phase[i]) +
                                       complex(0,1) * sin(phase[i])));

}

//==========================================================================

// Tau decay matrix element for tau decay into a single scalar meson.

// The maximum decay weight for this matrix element can be found analytically
// to be 4 * m_tau^2 * (m_tau^2 - m_meson^2). However, because m_tau >> m_meson
// for the relevant tau decay channels, this expression is approximated by
// m_tau^4.

//--------------------------------------------------------------------------

// Initialize constants for the helicity matrix element.

void HMETau2Meson::initConstants() {

  DECAYWEIGHTMAX = 4*pow4(pM[0]);

}

//--------------------------------------------------------------------------

// Initialize the hadronic current for the helicity matrix element.

void HMETau2Meson::initHadronicCurrent(vector<HelicityParticle>& p) {

  vector< Wave4 > u2;
  pMap[2] = 2;
  u2.push_back(Wave4(p[2].p()));
  u.push_back(u2);

}

//==========================================================================

// Tau decay matrix element for tau decay into two leptons. Because there is
// no hadronic current, but rather a leptonic current, the calculateME and
// initWaves methods must be redefined.

//--------------------------------------------------------------------------

// Initialize constants for the helicity matrix element.

void HMETau2TwoLeptons::initConstants() {

  DECAYWEIGHTMAX = 16*pow4(pM[0]);

}

//--------------------------------------------------------------------------

// Initialize spinors for the helicity matrix element.

void HMETau2TwoLeptons::initWaves(vector<HelicityParticle>& p) {

  u.clear();
  pMap.resize(4);
  setFermionLine(0,p[0],p[1]);
  setFermionLine(2,p[2],p[3]);

}

//--------------------------------------------------------------------------

// Return element for the helicity matrix element.

complex HMETau2TwoLeptons::calculateME(vector<int> h) {

  complex answer(0,0);
  for (int mu = 0; mu <= 3; mu++) {
    answer += (u[1][h[pMap[1]]] * gamma[mu] * (1 - gamma[5])
      * u[0][h[pMap[0]]]) * gamma[4](mu,mu) * (u[3][h[pMap[3]]]
      * gamma[mu] * (1 - gamma[5]) * u[2][h[pMap[2]]]);
  }
  return answer;

}

//==========================================================================

// Tau decay matrix element for tau decay into two mesons through an
// intermediate vector meson. This matrix element is used for pi^0 + pi^-
// decays (rho resonances), K^0 + K^- decays (rho resonances), and eta + K^-
// decays (K^* resonances). Note that for the rho resonances the pi^0 + pi^-
// running width dominates while for the K^* resonances the pi^- + K^0 running
// width dominates.

// vecM: on-shell masses for the vector resonances
// vecG: on-shell widths for the vector resonances
// vecP: phases used to calculate vector resonance weights
// vecA: amplitudes used to calculate vector resonance weights
// vecW: vector resonance weights

//--------------------------------------------------------------------------

// Initialize constants for the helicity matrix element.

void HMETau2TwoMesonsViaVector::initConstants() {

  // Clear the vectors from previous decays.
  vecM.clear(); vecG.clear(); vecP.clear(); vecA.clear(); vecW.clear();

  // Decay through K^* resonances (eta + K^- decay).
  if (abs(pID[2]) == 221) {
    DECAYWEIGHTMAX = 10;
    pM[2] = particleDataPtr->m0(211); pM[3] = particleDataPtr->m0(311);
    vecM.push_back(0.8921); vecM.push_back(1.700);
    vecG.push_back(0.0513); vecG.push_back(0.235);
    vecP.push_back(0);      vecP.push_back(M_PI);
    vecA.push_back(1);      vecA.push_back(0.038);
  }

  // Decay through rho resonances (pi^0 + pi^- and K^0 + K^- decays).
  else {
    if (abs(pID[2]) == 111)      DECAYWEIGHTMAX = 800;
    else if (abs(pID[2]) == 311) DECAYWEIGHTMAX = 6;
    pM[2] = particleDataPtr->m0(111); pM[3] = particleDataPtr->m0(211);
    vecM.push_back(0.7746); vecM.push_back(1.4080); vecM.push_back(1.700);
    vecG.push_back(0.1490); vecG.push_back(0.5020); vecG.push_back(0.235);
    vecP.push_back(0);      vecP.push_back(M_PI);   vecP.push_back(0);
    vecA.push_back(1.0);    vecA.push_back(0.167);  vecA.push_back(0.050);
  }
  calculateResonanceWeights(vecP, vecA, vecW);

}

//--------------------------------------------------------------------------

// Initialize the hadronic current for the helicity matrix element.

void HMETau2TwoMesonsViaVector::initHadronicCurrent(
  vector<HelicityParticle>& p) {

  vector< Wave4 > u2;
  Wave4 u3(p[3].p() - p[2].p());
  Wave4 u4(p[2].p() + p[3].p());
  double s1 = m2(u3, u4);
  double s2 = m2(u4);
  complex sumBW = 0;
  for (unsigned int i = 0; i < vecW.size(); i++)
    sumBW += vecW[i] * pBreitWigner(pM[2], pM[3], s2, vecM[i], vecG[i]);
  u2.push_back((u3 - s1 / s2 * u4) * sumBW);
  u.push_back(u2);

}

//==========================================================================

// Tau decay matrix element for tau decay into two mesons through both
// intermediate vector and scalar mesons.

// scaC: scalar coupling constant
// scaM: on-shell masses for the scalar resonances
// scaG: on-shell widths for the scalar resonances
// scaP: phases used to calculate scalar resonance weights
// scaA: amplitudes used to calculate scalar resonance weights
// scaW: scalar resonance weights
// vecC: scalar coupling constant
// vecM: on-shell masses for the vector resonances
// vecG: on-shell widths for the vector resonances
// vecP: phases used to calculate vector resonance weights
// vecA: amplitudes used to calculate vector resonance weights
// vecW: vector resonance weights

//--------------------------------------------------------------------------

// Initialize constants for the helicity matrix element.

void HMETau2TwoMesonsViaVectorScalar::initConstants() {

  DECAYWEIGHTMAX = 5400;
  // Clear the vectors from previous decays.
  scaM.clear(); scaG.clear(); scaP.clear(); scaA.clear(); scaW.clear();
  vecM.clear(); vecG.clear(); vecP.clear(); vecA.clear(); vecW.clear();
  // Scalar resonance parameters.
  scaC = 0.465;
  scaM.push_back(0.878);
  scaG.push_back(0.499);
  scaP.push_back(0);
  scaA.push_back(1);
  calculateResonanceWeights(scaP, scaA, scaW);
  // Vector resonance parameters.
  vecC = 1;
  vecM.push_back(0.89547); vecM.push_back(1.414);
  vecG.push_back(0.04619); vecG.push_back(0.232);
  vecP.push_back(0);       vecP.push_back(1.4399);
  vecA.push_back(1);       vecA.push_back(0.075);
  calculateResonanceWeights(vecP, vecA, vecW);

}

//--------------------------------------------------------------------------

// Initialize the hadronic current for the helicity matrix element.

void HMETau2TwoMesonsViaVectorScalar::initHadronicCurrent(
  vector<HelicityParticle>& p) {

  vector< Wave4 > u2;
  Wave4 u3(p[3].p() - p[2].p());
  Wave4 u4(p[2].p() + p[3].p());
  double s1 = m2(u3,u4);
  double s2 = m2(u4);
  complex scaSumBW = 0; complex scaSumW = 0;
  complex vecSumBW = 0; complex vecSumW = 0; complex vecSumBWM = 0;
  for (unsigned int i = 0; i < scaW.size(); i++) {
    scaSumBW  += scaW[i] * sBreitWigner(pM[2], pM[3], s2, scaM[i], scaG[i]);
    scaSumW   += scaW[i];
  }
  for (unsigned int i = 0; i < vecW.size(); i++) {
    vecSumBW  += vecW[i] * pBreitWigner(pM[2], pM[3], s2, vecM[i], vecG[i]);
    vecSumBWM += vecW[i] * pBreitWigner(pM[2], pM[3], s2, vecM[i], vecG[i]) /
        pow2(vecM[i]);
    vecSumW   += vecW[i];
  }
  u2.push_back(vecC * (vecSumBW * u3 - s1 * vecSumBWM * u4) / vecSumW +
                 scaC * u4 * scaSumBW / scaSumW);
  u.push_back(u2);

}

//==========================================================================

// Tau decay matrix element for tau decay into three mesons. This matrix
// element provides a base class for all implemented three meson decays.

// mode: three meson decay mode of the tau
// initMode(): initialize the decay mode
// initResonances(): initialize the resonance constants
// s1, s2, s3, s4: center-of-mass energies
// q, q2, q3, q4: summed and individual hadronic momentum four-vectors
// a1BW: stored value of a1BreitWigner for speed
// a1PhaseSpace(s): phase space factor for the a1
// a1BreitWigner(s): Breit-Wigner for the a1
// T(m0, m1, s, M, G, W): sum weighted p-wave Breit-Wigners
// T(s, M, G, W): sum weighted fixed width Breit-Wigners
// F1(), F2(), F3(), F4(): sub-current form factors

//--------------------------------------------------------------------------

// Initialize constants for the helicity matrix element.

void HMETau2ThreeMesons::initConstants() {

  initMode();
  initResonances();

}

//--------------------------------------------------------------------------

// Initialize the hadronic current for the helicity matrix element.

void HMETau2ThreeMesons::initHadronicCurrent(vector<HelicityParticle>& p) {

  vector< Wave4 > u2;

  // Initialize the momenta.
  initMomenta(p);

  // Calculate the center of mass energies.
  s1 = m2(q);
  s2 = m2(q3 + q4);
  s3 = m2(q2 + q4);
  s4 = m2(q2 + q3);

  // Calculate the form factors.
  a1BW = a1BreitWigner(s1);
  complex f1   = F1();
  complex f2   = F2();
  complex f3   = F3();
  complex f4   = F4();

  // Calculate the hadronic current.
  Wave4 u3 = (f3 - f2) * q2 + (f1 - f3) * q3 + (f2 - f1) * q4;
  u3 = u3 - (u3 * gamma[4] * q / s1) * q;
  if (f4 != complex(0, 0))
    u3 = u3 + complex(0, 1) * f4 * epsilon(q2, q3, q4);
  u2.push_back(u3);
  u.push_back(u2);

}

//--------------------------------------------------------------------------

// Initialize the tau decay mode.

void HMETau2ThreeMesons::initMode() {

  if (abs(pID[2]) == 111 && abs(pID[3]) == 111 && abs(pID[4]) == 211)
    mode = Pi0Pi0Pim;
  else if (abs(pID[2]) == 211 && abs(pID[3]) == 211 && abs(pID[4]) == 211)
    mode = PimPimPip;
  else if (abs(pID[2]) == 111 && abs(pID[3]) == 211 && abs(pID[4]) == 311)
    mode = Pi0PimK0b;
  else if (abs(pID[2]) == 211 && abs(pID[3]) == 211 && abs(pID[4]) == 321)
    mode = PimPipKm;
  else if (abs(pID[2]) == 111 && abs(pID[3]) == 211 && abs(pID[4]) == 221)
    mode = Pi0PimEta;
  else if (abs(pID[2]) == 211 && abs(pID[3]) == 321 && abs(pID[4]) == 321)
    mode = PimKmKp;
  else if (abs(pID[2]) == 111 && abs(pID[3]) == 311 && abs(pID[4]) == 321)
    mode = Pi0K0Km;
  else if (abs(pID[2]) == 130 && abs(pID[3]) == 211 && abs(pID[4]) == 310)
    mode = KlPimKs;
  else if (abs(pID[2]) == 111 && abs(pID[3]) == 111 && abs(pID[4]) == 321)
    mode = Pi0Pi0Km;
  else if (abs(pID[2]) == 130 && abs(pID[3]) == 130 && abs(pID[4]) == 211)
    mode = KlKlPim;
  else if (abs(pID[2]) == 211 && abs(pID[3]) == 310 && abs(pID[4]) == 310)
    mode = PimKsKs;
  else if (abs(pID[2]) == 211 && abs(pID[3]) == 311 && abs(pID[4]) == 311)
    mode = PimK0bK0;
  else
    mode = Uknown;
}

//--------------------------------------------------------------------------

// Initialize the momenta for the helicity matrix element.

void HMETau2ThreeMesons::initMomenta(vector<HelicityParticle>& p) {

  q = Wave4(p[2].p() + p[3].p() + p[4].p());
  // pi-, pi-, pi+ decay and pi0, pi0, pi- decay.
  if (mode == PimPimPip || mode == Pi0Pi0Pim) {
    q2 = Wave4(p[2].p()); q3 = Wave4(p[3].p()); q4 = Wave4(p[4].p());
  // K-, pi-, K+ decay.
  } else if (mode == PimKmKp) {
    q2 = Wave4(p[3].p()); q3 = Wave4(p[2].p()); q4 = Wave4(p[4].p());
  // K0, pi-, Kbar0 decay.
  } else if (mode == PimK0bK0) {
    q2 = Wave4(p[3].p()); q3 = Wave4(p[2].p()); q4 = Wave4(p[4].p());
  // K_S0, pi-, K_S0 decay.
  } else if (mode == PimKsKs) {
    q2 = Wave4(p[3].p()); q3 = Wave4(p[2].p()); q4 = Wave4(p[4].p());
  // K_L0, pi-, K_L0 decay.
  } else if (mode == KlKlPim) {
    q2 = Wave4(p[2].p()); q3 = Wave4(p[4].p()); q4 = Wave4(p[3].p());
  // K_S0, pi-, K_L0 decay.
  } else if (mode == KlPimKs) {
    q2 = Wave4(p[4].p()); q3 = Wave4(p[3].p()); q4 = Wave4(p[2].p());
  } // K-, pi0, K0 decay.
  else if (mode == Pi0K0Km) {
    q2 = Wave4(p[4].p()); q3 = Wave4(p[2].p()); q4 = Wave4(p[3].p());
  } // pi0, pi0, K- decay.
  else if (mode == Pi0Pi0Km) {
    q2 = Wave4(p[2].p()); q3 = Wave4(p[3].p()); q4 = Wave4(p[4].p());
  } // K-, pi-, pi+ decay.
  else if (mode == PimPipKm) {
    q2 = Wave4(p[4].p()); q3 = Wave4(p[2].p()); q4 = Wave4(p[3].p());
  } // pi-, Kbar0, pi0 decay.
  else if (mode == Pi0PimK0b) {
    q2 = Wave4(p[3].p()); q3 = Wave4(p[4].p()); q4 = Wave4(p[2].p());
  } // pi-, pi0, eta decay.
  else if (mode == Pi0PimEta) {
    q2 = Wave4(p[3].p()); q3 = Wave4(p[2].p()); q4 = Wave4(p[4].p());
  }

}

//--------------------------------------------------------------------------

// Return the phase space factor for the a1. Implements equation 3.16 of Z.
// Phys. C48 (1990) 445-452.

double HMETau2ThreeMesons::a1PhaseSpace(double s) {

  double piM  = 0.13957; // Mass of the charged pion.
  double rhoM = 0.773;   // Mass of the rho.
  if (s < pow2(3 * piM))
    return 0;
  else if (s < pow2(rhoM + piM)) {
    double sum = (s - 9 * piM * piM);
    return 4.1 * sum * sum * sum * (1 - 3.3 * sum + 5.8 * sum * sum);
  }
  else
    return s * (1.623 + 10.38 / s - 9.32 / (s * s) + 0.65 / (s * s * s));

}

//--------------------------------------------------------------------------

// Return the Breit-Wigner for the a1. Implements equation 3.18
// of Z. Phys. C48 (1990) 445-452.

complex HMETau2ThreeMesons::a1BreitWigner(double s) {

  double a1M = 1.251; // Mass of the a1.
  double a1G = 0.475; // Width of the a1.
  return a1M * a1M / (a1M * a1M - s - complex(0,1) * a1M * a1G
                      * a1PhaseSpace(s) / a1PhaseSpace(a1M * a1M));

}

//--------------------------------------------------------------------------

// Return summed weighted running p Breit-Wigner resonances.

complex HMETau2ThreeMesons::T(double m0, double m1, double s,
  vector<double> &M, vector<double> &G, vector<double> &W) {

  complex num(0, 0);
  double  den(0);
  for (unsigned int i = 0; i < M.size(); i++) {
    num += W[i] * pBreitWigner(m0, m1, s, M[i], G[i]);
    den += W[i];
  }
  return num / den;

}

//--------------------------------------------------------------------------

// Return summed weighted fixed width Breit-Wigner resonances.

complex HMETau2ThreeMesons::T(double s, vector<double> &M,
  vector<double> &G, vector<double> &W) {

  complex num(0, 0);
  double  den(0);
  for (unsigned int i = 0; i < M.size(); i++) {
    num += W[i] * breitWigner(s, M[i], G[i]);
    den += W[i];
  }
  return num / den;

}

//==========================================================================

// Tau decay matrix element for tau decay into three pions. This matrix element
// is taken from the Herwig++ implementation based on the CLEO fits.

// rhoM: on-shell masses for the rho resonances
// rhoG: on-shell widths for the rho resonances
// rhoPp: p-wave phase for the rho coupling to the a1
// rhoAp: p-wave amplitude for the rho coupling to the a1
// rhoPd: d-wave phase for the rho coupling to the a1
// rhoAd: d-wave amplitude for the rho coupling to the a1
// f0M: f0 on-shell mass
// f0G: f0 on-shell width
// f0P: phase for the coupling of the f0 to the a1
// f0A: amplitude for the coupling of the f0 to the a1
// f2M: f2 on-shell mass
// f2G: f2 on-shell width
// f2P: phase for the coupling of the f2 to the a1
// f2P: amplitude for the coupling of the f2 to the a1
// sigM: sigma on-shell mass
// sigG: sigma on-shell width
// sigP: phase for the coupling of the sigma to the a1
// sigA: amplitude for the coupling of the sigma to the a1

//--------------------------------------------------------------------------

// Initialize resonance constants for the helicity matrix element.

void HMETau2ThreePions::initResonances() {

  // Three charged pion decay.
  if (mode == PimPimPip) DECAYWEIGHTMAX = 6000;

  // Two neutral and one charged pion decay.
  else DECAYWEIGHTMAX = 3000;

  // Clear the vectors from previous decays.
  rhoM.clear(); rhoG.clear();
  rhoPp.clear(); rhoAp.clear(); rhoWp.clear();
  rhoPd.clear(); rhoAd.clear(); rhoWd.clear();

  // Rho parameters.
  rhoM.push_back(.7743);      rhoM.push_back(1.370);    rhoM.push_back(1.720);
  rhoG.push_back(.1491);      rhoG.push_back(.386);     rhoG.push_back(.250);
  rhoPp.push_back(0);         rhoPp.push_back(3.11018); rhoPp.push_back(0);
  rhoAp.push_back(1);         rhoAp.push_back(0.12);    rhoAp.push_back(0);
  rhoPd.push_back(-0.471239); rhoPd.push_back(1.66504); rhoPd.push_back(0);
  rhoAd.push_back(0.37);      rhoAd.push_back(0.87);    rhoAd.push_back(0);

  // Scalar and tensor parameters.
  f0M  = 1.186;    f2M  = 1.275;   sigM = 0.860;
  f0G  = 0.350;    f2G  = 0.185;   sigG = 0.880;
  f0P  = -1.69646; f2P  = 1.75929; sigP = 0.722566;
  f0A  = 0.77;     f2A  = 0.71;    sigA = 2.1;

  // Calculate the weights from the phases and amplitudes.
  calculateResonanceWeights(rhoPp, rhoAp, rhoWp);
  calculateResonanceWeights(rhoPd, rhoAd, rhoWd);
  f0W  = f0A  * (cos(f0P)  + complex(0,1) * sin(f0P));
  f2W  = f2A  * (cos(f2P)  + complex(0,1) * sin(f2P));
  sigW = sigA * (cos(sigP) + complex(0,1) * sin(sigP));

}

//--------------------------------------------------------------------------

// Return the first form factor.

complex HMETau2ThreePions::F1() {

  complex answer(0,0);

  // Three charged pion decay.
  if (mode == PimPimPip) {
    for (unsigned int i = 0; i < rhoM.size(); i++) {
      answer += - rhoWp[i] * pBreitWigner(pM[3], pM[4], s2, rhoM[i], rhoG[i])
        - rhoWd[i] / 3.0 * pBreitWigner(pM[2], pM[4], s3, rhoM[i], rhoG[i])
        * (s2 - s4);
    }
    answer += -2.0 / 3.0 * (sigW * sBreitWigner(pM[2], pM[4], s3, sigM, sigG)
            + f0W * sBreitWigner(pM[2], pM[4], s3, f0M, f0G));
    answer += f2W * (0.5 * (s4 - s3)
            * dBreitWigner(pM[3], pM[4], s2, f2M, f2G)
            - 1.0 / (18 * s3) * (4 * pow2(pM[2]) - s3)
            * (s1 + s3 - pow2(pM[2]))
            * dBreitWigner(pM[2], pM[4], s3, f2M, f2G));
  }

  // Two neutral and one charged pion decay.
  else {
    for (unsigned int i = 0; i < rhoM.size(); i++) {
      answer += rhoWp[i] * pBreitWigner(pM[3], pM[4], s2, rhoM[i], rhoG[i])
        - rhoWd[i] / 3.0 * pBreitWigner(pM[2], pM[4], s3, rhoM[i], rhoG[i])
        * (s4 - s2 - pow2(pM[4]) + pow2(pM[2]));
    }
    answer += 2.0 / 3.0 * (sigW * sBreitWigner(pM[2], pM[3], s4, sigM, sigG)
      + f0W * sBreitWigner(pM[2], pM[3], s4, f0M, f0G));
    answer += f2W / (18 * s4) * (s1 - pow2(pM[4]) + s4)
      * (4 * pow2(pM[2]) - s4) * dBreitWigner(pM[2], pM[3], s4, f2M, f2G);
  }
  return a1BW * answer;

}

//--------------------------------------------------------------------------

// Return the second form factor.

complex HMETau2ThreePions::F2() {

  complex answer(0,0);

  // Three charged pion decay.
  if (mode == PimPimPip) {
    for (unsigned int i = 0; i  < rhoM.size(); i++) {
      answer += -rhoWp[i] * pBreitWigner(pM[2], pM[4], s3, rhoM[i], rhoG[i])
        - rhoWd[i] / 3.0 * pBreitWigner(pM[3], pM[4], s2, rhoM[i], rhoG[i])
        * (s3 - s4);
    }
    answer += -2.0 / 3.0 * (sigW * sBreitWigner(pM[3], pM[4], s2, sigM, sigG)
      + f0W * sBreitWigner(pM[3], pM[4], s2, f0M, f0G));
    answer += f2W * (0.5 * (s4 - s2)
      * dBreitWigner(pM[2], pM[4], s3, f2M, f2G)
      - 1.0 / (18 * s2) * (4 * pow2(pM[2]) - s2) * (s1 + s2 - pow2(pM[2]))
      * dBreitWigner(pM[3], pM[4], s2, f2M, f2G));
  }

  // Two neutral and one charged pion decay.
  else {
    for (unsigned int i = 0; i < rhoM.size(); i++) {
        answer += -rhoWp[i] / 3.0 *
          pBreitWigner(pM[2], pM[4], s3, rhoM[i], rhoG[i]) -
          rhoWd[i] * pBreitWigner(pM[3], pM[4], s2, rhoM[i], rhoG[i]) *
          (s4 - s3 - pow2(pM[4]) + pow2(pM[3]));
    }
    answer += 2.0 / 3.0 * (sigW * sBreitWigner(pM[2], pM[3], s4, sigM, sigG)
                             + f0W * sBreitWigner(pM[2], pM[3], s4, f0M, f0G));
    answer += f2W / (18 * s4) * (s1 - pow2(pM[4]) + s4) *
        (4 * pow2(pM[2]) - s4) * dBreitWigner(pM[2], pM[3], s4, f2M, f2G);
  }
  return -a1BW * answer;

}

//--------------------------------------------------------------------------

// Return the third form factor.

complex HMETau2ThreePions::F3() {

  complex answer(0,0);

  // Three charged pion decay.
  if (mode == PimPimPip) {
    for (unsigned int i = 0; i < rhoM.size(); i++) {
      answer += -rhoWd[i] * (1.0 / 3.0 * (s3 - s4)
        * pBreitWigner(pM[3], pM[4], s2, rhoM[i], rhoG[i]) - 1.0 / 3.0
        * (s2 - s4) * pBreitWigner(pM[2], pM[4], s3, rhoM[i], rhoG[i]));
    }
    answer += -2.0 / 3.0 * (sigW * sBreitWigner(pM[3], pM[4], s2, sigM, sigG)
      + f0W * sBreitWigner(pM[3], pM[4], s2, f0M, f0G));
    answer += 2.0 / 3.0 * (sigW * sBreitWigner(pM[2], pM[4], s3, sigM, sigG)
      + f0W * sBreitWigner(pM[2], pM[4], s3, f0M, f0G));
    answer += f2W * (-1.0 / (18 * s2) * (4 * pow2(pM[2]) - s2)
      * (s1 + s2 - pow2(pM[2])) * dBreitWigner(pM[3], pM[4], s2, f2M, f2G)
      + 1.0 / (18 * s3) * (4 * pow2(pM[2]) - s3) * (s1 + s3 - pow2(pM[2]))
      * dBreitWigner(pM[2], pM[4], s3, f2M, f2G));
  }

  // Two neutral and one charged pion decay.
  else {
    for (unsigned int i = 0; i < rhoM.size(); i++) {
      answer += rhoWd[i] * (-1.0 / 3.0 * (s4 - s3 - pow2(pM[4]) + pow2(pM[3]))
        * pBreitWigner(pM[3], pM[4], s2, rhoM[i], rhoG[i])
        + 1.0 / 3.0 * (s4 - s2 - pow2(pM[4]) + pow2(pM[2]))
        * pBreitWigner(pM[2], pM[4], s3, rhoM[i], rhoG[i]));
    }
    answer += -f2W / 2.0 * (s2 - s3)
      * dBreitWigner(pM[2], pM[3], s4, f2M, f2G);
  }
  return a1BW * answer;

}

//--------------------------------------------------------------------------

// Return the running width for the a1 (multiplied by a factor of a1M).

double HMETau2ThreePions::a1PhaseSpace(double s) {

  double picM = 0.1753; // (m_pi^- + m_pi^- + m_pi^+)^2
  double pinM = 0.1676; // (m_pi^0 + m_pi^0 + m_pi^-)^2
  double kM   = 0.496;  // K mass.
  double ksM  = 0.894;  // K^* mass.
  double picG = 0;      // Width contribution from three charged pions.
  double pinG = 0;      // Width contribution from two neutral one charged.
  double kG   = 0;      // Width contributions from s-wave K K^*.
  double piW  = pow2(0.2384)/1.0252088; // Overall weight factor.
  double kW   = pow2(4.7621);           // K K^* width weight factor.

  // Three charged pion width contribution.
  if (s < picM)
    picG = 0;
  else if (s < 0.823)
    picG = 5.80900 * pow3(s - picM) * (1.0 - 3.00980 * (s - picM) +
                                         4.5792 * pow2(s - picM));
  else
    picG = -13.91400 + 27.67900 * s - 13.39300 * pow2(s) + 3.19240 * pow3(s)
        - 0.10487 * pow4(s);

  // Two neutral and one charged pion width contribution.
  if (s < pinM)
    pinG = 0;
  else if (s < 0.823)
    pinG = 6.28450 * pow3(s - pinM) * (1.0 - 2.95950 * (s - pinM) +
                                         4.33550 * pow2(s - pinM));
  else
    pinG = -15.41100 + 32.08800 * s - 17.66600 * pow2(s) + 4.93550 * pow3(s)
        - 0.37498 * pow4(s);

  // K and K^* width contribution.
  if (s > pow2(ksM + kM))
    kG = 0.5 * sqrt((s - pow2(ksM + kM)) * (s - pow2(ksM - kM))) / s;
  return piW*(picG + pinG + kW*kG);

}

//--------------------------------------------------------------------------

// Return the Breit-Wigner for the a1.

complex HMETau2ThreePions::a1BreitWigner(double s) {

  double a1M = 1.331; // Mass of the a1.
  return a1M*a1M/(a1M*a1M - s - complex(0,1)*a1PhaseSpace(s));

}

//==========================================================================

// Tau decay matrix element for tau decay into three mesons with kaons.
// The form factors are taken from hep-ph/9503474.

// rhoMa(v): on-shell masses for the axial (vector) rho resonances
// rhoGa(v): widths for the axial (vector) rho resonances
// rhoWa(v): weights for the axial (vector) rho resonances
// kstarXa(v): on-shell masses, widths, and weights for the K* resonances
// k1Xa(b): on-shell masses, width, and weight for the K1 resonances,
//          the a constants are for K1 -> K* pi, K* -> K pi
//          the b constants are for K1 -> rho K, rho -> pi pi
// omegaX: on-shell masses, width, and weights for the omega reosnances
// kM:  kaon mass
// piM: charged pion mass
// piW: pion coupling, f_W

//--------------------------------------------------------------------------

// Initialize resonance constants for the helicity matrix element.

void HMETau2ThreeMesonsWithKaons::initResonances() {

  // K-, pi-, K+ decay.
  if (mode == PimKmKp) DECAYWEIGHTMAX = 130;
  // K0, pi-, Kbar0 decay.
  else if (mode == PimK0bK0) DECAYWEIGHTMAX = 115;
  // K_S0, pi-, K_S0 and K_L0, pi-, K_L0 decay.
  else if (mode == PimKsKs || mode == KlKlPim) DECAYWEIGHTMAX = 230;
  // K_S0, pi-, K_L0 decay.
  else if (mode == KlPimKs) DECAYWEIGHTMAX = 230;
  // K-, pi0, K0 decay.
  else if (mode == Pi0K0Km) DECAYWEIGHTMAX = 125;
  // pi0, pi0, K- decay.
  else if (mode == Pi0Pi0Km) DECAYWEIGHTMAX = 2.5e4;
  // K-, pi-, pi+ decay.
  else if (mode == PimPipKm) DECAYWEIGHTMAX = 1.8e4;
  // pi-, Kbar0, pi0 decay.
  else if (mode == Pi0PimK0b) DECAYWEIGHTMAX = 3.9e4;

  // Clear the vectors from previous decays.
  rhoMa.clear();   rhoGa.clear();   rhoWa.clear();
  rhoMv.clear();   rhoGv.clear();   rhoWv.clear();
  kstarMa.clear(); kstarGa.clear(); kstarWa.clear();
  kstarMv.clear(); kstarGv.clear(); kstarWv.clear();
  k1Ma.clear();    k1Ga.clear();    k1Wa.clear();
  k1Mb.clear();    k1Gb.clear();    k1Wb.clear();
  omegaM.clear();  omegaG.clear();  omegaW.clear();

  // Rho parameters.
  rhoMa.push_back(0.773); rhoGa.push_back(0.145); rhoWa.push_back(1);
  rhoMa.push_back(1.370); rhoGa.push_back(0.510); rhoWa.push_back(-0.145);
  rhoMv.push_back(0.773); rhoGv.push_back(0.145); rhoWv.push_back(1);
  rhoMv.push_back(1.500); rhoGv.push_back(0.220); rhoWv.push_back(-6.5 / 26.0);
  rhoMv.push_back(1.750); rhoGv.push_back(0.120); rhoWv.push_back(-1.0 / 26.0);

  // Kstar parameters.
  kstarMa.push_back(0.892); kstarGa.push_back(0.050);
  kstarMa.push_back(1.412); kstarGa.push_back(0.227);
  kstarWa.push_back(1);
  kstarWa.push_back(-0.135);
  kstarMv.push_back(0.892); kstarGv.push_back(0.050);
  kstarMv.push_back(1.412); kstarGv.push_back(0.227);
  kstarMv.push_back(1.714); kstarGv.push_back(0.323);
  kstarWv.push_back(1);
  kstarWv.push_back(-6.5 / 26.0);
  kstarWv.push_back(-1.0 / 26.0);

  // K1 parameters.
  k1Ma.push_back(1.270); k1Ga.push_back(0.090); k1Wa.push_back(0.33);
  k1Ma.push_back(1.402); k1Ga.push_back(0.174); k1Wa.push_back(1);
  k1Mb.push_back(1.270); k1Gb.push_back(0.090); k1Wb.push_back(1);

  // Omega and phi parameters.
  omegaM.push_back(0.782); omegaG.push_back(0.00843); omegaW.push_back(1);
  omegaM.push_back(1.020); omegaG.push_back(0.00443); omegaW.push_back(0.05);

  // Kaon and pion parameters
  kM = 0.49765; piM = 0.13957; piW = 0.0942;

}

//--------------------------------------------------------------------------

// Return the first form factor.

complex HMETau2ThreeMesonsWithKaons::F1() {

  complex answer;
  // K-, pi-, K+ decay.
  if (mode == PimKmKp)
    answer = a1BW * T(piM, kM, s2, kstarMa, kstarGa, kstarWa) / 2.0;
  // K0, pi-, Kbar0 decay.
  else if (mode == PimK0bK0)
    answer = a1BW * T(piM, kM, s2, kstarMa, kstarGa, kstarWa) / 2.0;
  // K_S0, pi-, K_S0 decay and K_L0, pi-, K_L0 decay.
  else if (mode == PimKsKs || mode == KlKlPim)
    answer = -a1BW * (T(piM, kM, s2, kstarMa, kstarGa, kstarWa)
                      + T(piM, kM, s4, kstarMa, kstarGa, kstarWa)) / 2.0;
  // K_S0, pi-, K_L0 decay.
  else if (mode == KlPimKs)
    answer = a1BW * (T(piM, kM, s2, kstarMa, kstarGa, kstarWa)
                      - T(piM, kM, s4, kstarMa, kstarGa, kstarWa)) / 2.0;
  // K-, pi0, K0 decay.
  else if (mode == Pi0K0Km)
    answer = a1BW * (T(piM, kM, s2, kstarMa, kstarGa, kstarWa)
                     - T(piM, kM, s4, kstarMa, kstarGa, kstarWa)) / 2.0;
  // pi0, pi0, K- decay.
  else if (mode == Pi0Pi0Km)
    answer = T(s1, k1Ma, k1Ga, k1Wa)
      * T(piM, kM, s2, kstarMa, kstarGa, kstarWa);
  // K-, pi-, pi+ decay.
  else if (mode == PimPipKm)
    answer = T(s1, k1Mb, k1Gb, k1Wb)
      * T(piM, piM, s2, rhoMa, rhoGa, rhoWa);
  // pi-, Kbar0, pi0 decay.
  else if (mode == Pi0PimK0b)
    answer = T(s1, k1Ma, k1Ga, k1Wa)
      * (T(piM, kM, s2, kstarMa, kstarGa, kstarWa)
         - T(piM, kM, s4, kstarMa, kstarGa, kstarWa));
  return -1.0 / 3.0 * answer;
}

//--------------------------------------------------------------------------

// Return the second form factor.

complex HMETau2ThreeMesonsWithKaons::F2() {

  complex answer;
  // K-, pi-, K+ decay.
  if (mode == PimKmKp)
    answer =  a1BW * T(piM, piM, s3, rhoMa, rhoGa, rhoWa) / 2.0;
  // K0, pi-, Kbar0 decay.
  else if (mode == PimK0bK0)
    answer =  a1BW * T(piM, piM, s3, rhoMa, rhoGa, rhoWa) / 2.0;
  // K_S0, pi-, K_S0 decay and K_L0, pi-, K_L0 decay.
  else if (mode == PimKsKs || mode == KlKlPim)
    answer = a1BW * T(piM, kM, s4, kstarMa, kstarGa, kstarWa) / 2.0;
  // K_S0, pi-, K_L0 decay.
  else if (mode == KlPimKs)
    answer = a1BW * (2.0 * T(piM, piM, s3, rhoMa, rhoGa, rhoWa)
                     + T(piM, kM, s4, kstarMa, kstarGa, kstarWa)) / 2.0;
  // K-, pi0, K0 decay.
  else if (mode == Pi0K0Km)
    answer = a1BW * (2.0 * T(piM, piM, s3, rhoMa, rhoGa, rhoWa)
                     + T(piM, kM, s4, kstarMa, kstarGa, kstarWa)) / 2.0;
  // pi0, pi0, K- decay.
  else if (mode == Pi0Pi0Km)
    answer = T(s1, k1Ma, k1Ga, k1Wa)
      * T(piM, kM, s3, kstarMa, kstarGa, kstarWa);
  // K-, pi-, pi+ decay.
  else if (mode == PimPipKm)
    answer = T(s1, k1Ma, k1Ga, k1Wa)
      * T(piM, kM, s3, kstarMa, kstarGa, kstarWa);
  // pi-, Kbar0, pi0 decay.
  else if (mode == Pi0PimK0b)
    answer = 2.0 * T(s1, k1Mb, k1Gb, k1Wb)
      * T(piM, piM, s3, rhoMa, rhoGa, rhoWa)
      + T(s1, k1Ma, k1Ga, k1Wa) * T(piM, kM, s4, kstarMa, kstarGa, kstarWa);
  return 1.0 / 3.0 * answer;

}

//--------------------------------------------------------------------------

// Return the fourth form factor.

complex HMETau2ThreeMesonsWithKaons::F4() {

  complex answer;
  // K-, pi-, K+ decay.
  if (mode == PimKmKp)
    answer = (sqrt(2.) - 1) * T(piM, piM, s1, rhoMv, rhoGv, rhoWv)
      * (sqrt(2.) * T(s3, omegaM, omegaG, omegaW)
         + T(piM, kM, s2, kstarMa, kstarGa, kstarWa));
  // K0, pi-, Kbar0 decay.
  else if (mode == PimK0bK0)
    answer = -(sqrt(2.) - 1) * T(piM, piM, s1, rhoMv, rhoGv, rhoWv)
      * (sqrt(2.) * T(s3, omegaM, omegaG, omegaW)
         + T(piM, kM, s2, kstarMa, kstarGa, kstarWa));
  // K_S0, pi-, K_S0 decay and K_L0, pi-, K_L0 decay.
  else if (mode == PimKsKs || mode == KlKlPim)
    answer = (sqrt(2.) - 1) * T(piM, piM, s1, rhoMv, rhoGv, rhoWv)
      * (T(piM, kM, s2, kstarMa, kstarGa, kstarWa)
         - T(piM, kM, s4, kstarMa, kstarGa, kstarWa));
  // K_S0, pi-, K_L0 decay.
  else if (mode == KlPimKs)
    answer = -(sqrt(2.) - 1) * T(piM, piM, s1, rhoMv, rhoGv, rhoWv)
      * (2 * sqrt(2.) * T(s3, omegaM, omegaG, omegaW)
         + T(piM, kM, s2, kstarMa, kstarGa, kstarWa)
         + T(piM, kM, s4, kstarMa, kstarGa, kstarWa));
  // K-, pi0, K0 decay.
  else if (mode == Pi0K0Km)
    answer = -(sqrt(2.) - 1) * T(piM, piM, s1, rhoMv, rhoGv, rhoWv)
      * (T(piM, kM, s4, kstarMa, kstarGa, kstarWa)
         - T(piM, kM, s2, kstarMa, kstarGa, kstarWa));
  // pi0, pi0, K- decay.
  else if (mode == Pi0Pi0Km)
    answer = T(piM, kM, s1, kstarMv, kstarGv, kstarWv)
      * (T(piM, kM, s2, kstarMa, kstarGa, kstarWa)
         - T(piM, kM, s3, kstarMa, kstarGa, kstarWa));
  // K-, pi-, pi+ decay.
  else if (mode == PimPipKm)
    answer = -T(piM, kM, s1, kstarMv, kstarGv, kstarWv)
      * (T(piM, piM, s2, rhoMa, rhoGa, rhoWa)
         + T(piM, kM, s3, kstarMa, kstarGa, kstarWa));
  // pi-, Kbar0, pi0 decay.
  else if (mode == Pi0PimK0b)
    answer = T(piM, kM, s1, kstarMv, kstarGv, kstarWv)
      * (2.0 * T(piM, piM, s3, rhoMa, rhoGa, rhoWa)
         + T(piM, kM, s2, kstarMa, kstarGa, kstarWa)
         + T(piM, kM, s4, kstarMa, kstarGa, kstarWa));
  return 1.0 / (8.0 * M_PI * M_PI * piW * piW) * answer;

}

//==========================================================================

// Tau decay matrix element for tau decay into three mesons. The form
// factors are taken from Comp. Phys. Com. 76 (1993) 361-380.

// rhoMa(v): on-shell masses for the axial (vector) rho resonances
// rhoGa(v): widths for the axial (vector) rho resonances
// rhoWa(v): weights for the axial (vector) rho resonances
// kstarX: on-shell masses, widths, and weights for the K* resonances
// k1X: on-shell masses, width, and weight for the K1 resonances
// kM:  kaon mass
// piM: charged pion mass
// piW: pion coupling, f_W

//--------------------------------------------------------------------------

// Initialize resonances for the helicity matrix element.

void HMETau2ThreeMesonsGeneric::initResonances() {

  // pi-, pi-, pi+ decay and pi0, pi0, pi- decay.
  if (mode == PimPimPip || mode == Pi0Pi0Pim) DECAYWEIGHTMAX = 1.3e4;
  // K-, pi-, K+ decay.
  else if (mode == PimKmKp) DECAYWEIGHTMAX = 330;
  // K0, pi-, Kbar0 decay.
  else if (mode == PimK0bK0) DECAYWEIGHTMAX = 300;
  // K-, pi0, K0 decay.
  else if (mode == Pi0K0Km) DECAYWEIGHTMAX = 40;
  // pi0, pi0, K- decay.
  else if (mode == Pi0Pi0Km) DECAYWEIGHTMAX = 9.4e4;
  // K-, pi-, pi+ decay.
  else if (mode == PimPipKm) DECAYWEIGHTMAX = 9.0e3;
  // pi-, Kbar0, pi0 decay.
  else if (mode == Pi0PimK0b) DECAYWEIGHTMAX = 1.2e4;
  // pi-, pi0, eta decay.
  else if (mode == Pi0PimEta) DECAYWEIGHTMAX = 360;

  // Clear the vectors from previous decays.
  rhoMa.clear();   rhoGa.clear();   rhoWa.clear();
  rhoMv.clear();   rhoGv.clear();   rhoWv.clear();
  kstarM.clear();  kstarG.clear();  kstarW.clear();
  k1M.clear();     k1G.clear();     k1W.clear();

  // Rho parameters.
  rhoMa.push_back(0.773); rhoGa.push_back(0.145); rhoWa.push_back(1);
  rhoMa.push_back(1.370); rhoGa.push_back(0.510); rhoWa.push_back(-0.145);
  rhoMv.push_back(0.773); rhoGv.push_back(0.145); rhoWv.push_back(-26);
  rhoMv.push_back(1.5);   rhoGv.push_back(0.220); rhoWv.push_back(6.5);
  rhoMv.push_back(1.75);  rhoGv.push_back(0.120); rhoWv.push_back(1);

  // Kaon parameters.
  kstarM.push_back(0.892); kstarG.push_back(0.0513); kstarW.push_back(1);
  k1M.push_back(1.402);    k1G.push_back(0.174);     k1W.push_back(1);

  // Kaon and pion parameters
  kM = 0.49765; piM = 0.13957; piW = 0.0942;

}

//--------------------------------------------------------------------------

// Return the first form factor.

complex HMETau2ThreeMesonsGeneric::F1() {

  complex answer;
  // pi-, pi-, pi+ decay and pi0, pi0, pi- decay.
  if (mode == PimPimPip || mode == Pi0Pi0Pim)
    answer = a1BW * T(piM, piM, s2, rhoMa, rhoGa, rhoWa);
  // K-, pi-, K+ decay.
  else if (mode == PimKmKp)
    answer = -a1BW * T(piM, kM, s2, kstarM, kstarG, kstarW) / 3.0;
  // K0, pi-, Kbar0 decay.
  else if (mode == PimK0bK0)
    answer = -a1BW * T(piM, kM, s2, kstarM, kstarG, kstarW) / 3.0;
  // K-, pi0, K0 decay.
  else if (mode == Pi0K0Km)
    answer = 0;
  // pi0, pi0, K- decay.
  else if (mode == Pi0Pi0Km)
    answer = T(s1, k1M, k1G, k1W) * T(piM, kM, s2, kstarM, kstarG, kstarW);
  // K-, pi-, pi+ decay.
  else if (mode == PimPipKm)
    answer = -T(s1, k1M, k1G, k1W) * T(piM, piM, s2, rhoMa, rhoGa, rhoWa)
           / 3.0;
  // pi-, Kbar0, pi0 decay.
  else if (mode == Pi0PimK0b)
    answer = 0;
  // pi-, pi0, eta decay.
  else if (mode == Pi0PimEta)
    answer = 0;
  return answer;

}

//--------------------------------------------------------------------------

// Return the second form factor.

complex HMETau2ThreeMesonsGeneric::F2() {

  complex answer;
  // pi-, pi-, pi+ decay and pi0, pi0, pi- decay.
  if (mode == PimPimPip || mode == Pi0Pi0Pim)
    answer =  -a1BW * T(piM, piM, s3, rhoMa, rhoGa, rhoWa);
  // K-, pi-, K+ decay.
  else if (mode == PimKmKp)
    answer = a1BW * T(piM, piM, s3, rhoMa, rhoGa, rhoWa) / 3.0;
  // K0, pi-, Kbar0 decay.
  else if (mode == PimK0bK0)
    answer = a1BW * T(piM, piM, s3, rhoMa, rhoGa, rhoWa) / 3.0;
  // K-, pi0, K0 decay.
  else if (mode == Pi0K0Km)
    answer = a1BW * T(piM, piM, s3, rhoMa, rhoGa, rhoWa);
  // pi0, pi0, K- decay.
  else if (mode == Pi0Pi0Km)
    answer = -T(s1, k1M, k1G, k1W) * T(piM, kM, s3, kstarM, kstarG, kstarW);
  // K-, pi-, pi+ decay.
  else if (mode == PimPipKm)
    answer = T(s1, k1M, k1G, k1W)
      * T(piM, kM, s3, kstarM, kstarG, kstarW) / 3.0;
  // pi-, Kbar0, pi0 decay.
  else if (mode == Pi0PimK0b)
    answer = T(s1, k1M, k1G, k1W) * T(piM, piM, s3, rhoMa, rhoGa, rhoWa);
  // pi-, pi0, eta decay.
  else if (mode == Pi0PimEta)
    answer = 0;
  return answer;

}

//--------------------------------------------------------------------------

// Return the fourth form factor.

complex HMETau2ThreeMesonsGeneric::F4() {

  complex answer;
  // pi-, pi-, pi+ decay and pi0, pi0, pi- decay.
  if (mode == PimPimPip || mode == Pi0Pi0Pim)
    answer = 0;
  // K-, pi-, K+ decay.
  else if (mode == PimKmKp)
    answer = T(piM, piM, s1, rhoMv, rhoGv, rhoWv)
      * (T(piM, piM, s3, rhoMa, rhoGa, rhoWa)
         - 0.2 * T(piM, kM, s2, kstarM, kstarG, kstarW)) * (1.25);
  // K0, pi-, Kbar0 decay.
  else if (mode == PimK0bK0)
    answer = -T(piM, piM, s1, rhoMv, rhoGv, rhoWv)
      * (T(piM, piM, s3, rhoMa, rhoGa, rhoWa)
         - 0.2 * T(piM, kM, s2, kstarM, kstarG, kstarW)) * (1.25);
  // K-, pi0, K0 decay.
  else if (mode == Pi0K0Km)
    answer = 0;
  // pi0, pi0, K- decay.
  else if (mode == Pi0Pi0Km)
    answer = 0;
  // K-, pi-, pi+ decay.
  else if (mode == PimPipKm)
    answer = -T(piM, kM, s1, kstarM, kstarG, kstarW)
      * (T(piM, piM, s2, rhoMa, rhoGa, rhoWa)
         - 0.2 * T(piM, kM, s3, kstarM, kstarG, kstarW)) * (1.25);
  // pi-, Kbar0, pi0 decay.
  else if (mode == Pi0PimK0b)
    answer = 2.0 * T(piM, kM, s1, kstarM, kstarG, kstarW)
      * (T(piM, piM, s3, rhoMa, rhoGa, rhoWa)
         - 0.2 * T(piM, kM, s2, kstarM, kstarG, kstarW)) * (1.25);
  // pi-, pi0, eta decay.
  else if (mode == Pi0PimEta)
    answer = T(piM, piM, s1, rhoMv, rhoGv, rhoWv)
      * T(piM, piM, s4, rhoMa, rhoGa, rhoWa);
  return 1.0 / (4.0 * M_PI * M_PI * piW * piW) * answer;

}

//==========================================================================

// Tau decay matrix element for tau decay into two pions with a photons taken
// from Comp. Phys. Com. 76 (1993) 361-380. Because a photon is in the final
// state the matrix element is reimplented to handle the two spin states.

// F(s, M, G, W): form factor for resonance
// rhoM: on-shell mass of the rho resonances
// rhoG: width of the rho resonances
// rhoW: weight of the rho resonances
// omegaX: on-shell mass, width, and weight of the omega resonances
// piM: charged pion mass

//--------------------------------------------------------------------------

// Initialize constants for the helicity matrix element.

void HMETau2TwoPionsGamma::initConstants() {

  DECAYWEIGHTMAX = 4e4;

  // Clear the vectors from previous decays.
  rhoM.clear();   rhoG.clear();   rhoW.clear();
  omegaM.clear(); omegaG.clear(); omegaW.clear();

  // Set parameters.
  rhoM.push_back(0.773);   rhoG.push_back(0.145);    rhoW.push_back(1);
  rhoM.push_back(1.7);     rhoG.push_back(0.26);     rhoW.push_back(-0.1);
  omegaM.push_back(0.782); omegaG.push_back(0.0085); omegaW.push_back(1);
  piM = 0.13957;

}

//--------------------------------------------------------------------------

// Initialize wave functions for the helicity matrix element.
void HMETau2TwoPionsGamma::initWaves(vector<HelicityParticle>& p) {

  // Calculate the hadronic current.
  u.clear();
  pMap.resize(p.size());
  setFermionLine(0, p[0], p[1]);

  // Calculate the hadronic current.
  vector< Wave4 > u2;
  Wave4 q(p[2].p() + p[3].p() + p[4].p());
  Wave4 q2(p[2].p()), q3(p[3].p()), q4(p[4].p());
  double s1 = m2(q);
  double s2 = m2(q3 + q2);
  complex f = F(s1, rhoM, rhoG, rhoW) * F(0, rhoM, rhoG, rhoW)
    * F(s2, omegaM, omegaG, omegaW);
  double q4q2 = m2(q4, q2);
  double q4q3 = m2(q4, q3);
  double q3q2 = m2(q3, q2);
  for (int h = 0; h < 2; h++) {
    Wave4 e = p[2].wave(h);
    complex q4e = q4*gamma[4]*e;
    complex q3e = q3*gamma[4]*e;
    u2.push_back(f * (e * (piM*piM*q4q2 - q3q2*(q4q3 - q4q2))
                      - q3 * (q3e*q4q2 - q4e*q3q2)
                      + q2 * (q3e*q4q3 - q4e*(piM*piM + q3q2))));
  }
  u.push_back(u2);

}

//--------------------------------------------------------------------------

// Return element for the helicity matrix element.
complex HMETau2TwoPionsGamma::calculateME(vector<int> h) {

  complex answer(0,0);
  for (int mu = 0; mu <= 3; mu++) {
    answer +=
        (u[1][h[pMap[1]]] * gamma[mu] * (1 - gamma[5]) * u[0][h[pMap[0]]])
        * gamma[4](mu,mu) * u[2][h[2]](mu);
  }
  return answer;

}

//--------------------------------------------------------------------------

// Return the form factor.
complex HMETau2TwoPionsGamma::F(double s, vector<double> M, vector<double> G,
                                vector<double> W) {

  complex answer(0, 0);
  for (unsigned int i = 0; i < M.size(); i++)
    answer += W[i] / (M[i]*M[i] - s - complex(0, 1) * M[i] * G[i]);
  return answer;

}

//==========================================================================

// Tau decay matrix element for tau decay into pions using the Novosibirsk
// model of Comp. Phys. Com. 174 (2006) 818-835.

// G(i,s): G-factor for index 1, 2, or 3
// tX(q,q1,q2,q3,q4): combined resonance current
// a1D(s): a1 Breit-Wigner denominator
// rhoD(s): rho Breit-Wigner denominator
// sigD(s): sigma Breit-Wigner denominator
// omeD(s): omega Breit-Wigner denominator
// a1FormFactor(s): form factor for the a1 resonance
// rhoFormFactor1(2)(s): form factor for the rho resonance
// omeFormFactor(s): form factor for the omega resonance
// sigM: on-shell mass of the sigma resonance
// sigG: width of the sigma resonance
// sigA: amplitude of the sigma resonance
// sigP: phase of the sigma resonance
// sigW: weight of the sigma resonance (from amplitude and weight)
// omeX: mass, width, amplitude, phase, and weight of the omega resonance
// a1X: mass and width of the a1 resonance
// rhoX: mass and width of the rho resonance
// picM: charge pion mass
// pinM: neutral pion mass
// lambda2: a1 form factor cut-off

//--------------------------------------------------------------------------

// Initialize constants for the helicity matrix element.

void HMETau2FourPions::initConstants() {

  if (abs(pID[3]) == 111) DECAYWEIGHTMAX = 5e8;
  else DECAYWEIGHTMAX = 5e9;
  pinM  = particleDataPtr->m0(111);
  picM  = particleDataPtr->m0(211);
  sigM = 0.8;     omeM = 0.782;   a1M  = 1.23; rhoM = 0.7761;
  sigG = 0.8;     omeG = 0.00841; a1G  = 0.45; rhoG = 0.1445;
  sigP = 0.43585; omeP = 0.0;
  sigA = 1.39987; omeA = 1.0;
  sigW = sigA*(cos(sigP)+complex(0,1)*sin(sigP));
  omeW = omeA*(cos(omeP)+complex(0,1)*sin(omeP));
  lambda2 = 1.2;

}

//--------------------------------------------------------------------------

// Initialize the hadronic current for the helicity matrix element.

void HMETau2FourPions::initHadronicCurrent(vector<HelicityParticle>& p) {

  vector< Wave4 > u2;

  // Create pion momenta.
  Wave4 q(p[2].p() + p[3].p() + p[4].p()+ p[5].p());
  Wave4 q2(p[2].p()), q3(p[3].p()), q4(p[4].p()), q5(p[5].p());

  // Calculate the four pion system energy.
  double s = m2(q);

  // Create the hadronic current for the 3 neutral pion channel.
  if (abs(pID[3]) == 111)
    u2.push_back(G(1,s)*(t1(q,q3,q4,q5,q2) + t1(q,q3,q2,q5,q4) +
                         t1(q,q4,q3,q5,q2) + t1(q,q4,q2,q5,q3) +
                         t1(q,q2,q3,q5,q4) + t1(q,q2,q4,q5,q3) +
                         t2(q,q3,q5,q4,q2) + t2(q,q4,q5,q3,q2) +
                         t2(q,q2,q5,q4,q3) - t2(q,q5,q3,q4,q2) -
                         t2(q,q5,q4,q3,q2) - t2(q,q5,q2,q4,q3)));

  // Create the hadronic current for the 3 charged pion channel.
  else if (abs(pID[3]) == 211)
    u2.push_back(G(2,s)*(t1(q,q3,q5,q4,q2) + t1(q,q4,q5,q3,q2) +
                         t1(q,q3,q4,q5,q2) + t1(q,q4,q3,q5,q2) +
                         t1(q,q2,q4,q3,q5) + t1(q,q2,q3,q4,q5) +
                         t2(q,q2,q4,q3,q5) + t2(q,q2,q3,q4,q5) -
                         t2(q,q3,q2,q4,q5) - t2(q,q4,q2,q3,q5)) +
                 G(3,s)*(t3(q,q3,q5,q4,q2) + t3(q,q4,q5,q3,q2) -
                         t3(q,q3,q4,q5,q2) - t3(q,q4,q3,q5,q2) -
                         t3(q,q3,q2,q4,q5) - t3(q,q4,q2,q3,q5)));
  u.push_back(u2);

}

//--------------------------------------------------------------------------

// Return the first t-vector.

Wave4 HMETau2FourPions::t1(Wave4 &q, Wave4 &q1, Wave4 &q2,
                           Wave4 &q3, Wave4 &q4) {

  Wave4  a1Q(q2 + q3 + q4);
  Wave4 rhoQ(q3 + q4);
  double  a1S = m2(a1Q);
  double rhoS = m2(rhoQ);

  // Needed to match Herwig++.
  double gM = sqrtpos(rhoM*rhoM - 4*picM*picM) * (rhoM*rhoM - 4*picM*picM)
    / rhoM;
  double dm = (rhoFormFactor1(0) - rhoFormFactor1(rhoM*rhoM) +
                 rhoM*rhoM * rhoFormFactor2(rhoM*rhoM)) / gM;
  return - a1FormFactor(a1S) / (a1D(a1S) * rhoD(rhoS)) * pow2(a1M) *
    (rhoM*rhoM + rhoM*rhoG*dm) *
    (m2(q,a1Q) *  (m2(q3,a1Q) * q4 - m2(q4,a1Q) * q3) +
     (m2(q,q4) * m2(q1,q3) - m2(q,q3) * m2(q1,q4)) * a1Q);

}

//--------------------------------------------------------------------------

// Return the second t-vector.

Wave4 HMETau2FourPions::t2(Wave4 &q, Wave4 &/*q1*/, Wave4 &q2,
                           Wave4 &q3, Wave4 &q4) {

  Wave4  a1Q(q2 + q3 + q4);
  Wave4 sigQ(q3 + q4);
  double  a1S = m2(a1Q);
  double sigS = m2(sigQ);
  return sigW * a1FormFactor(a1S) / (a1D(a1S) * sigD(sigS)) *
    pow2(a1M) * pow2(sigM) * (m2(q,a1Q) * a1S * q2 - m2(q,q2) * a1S * a1Q);

}

//--------------------------------------------------------------------------

// Return the third t-vector.

Wave4 HMETau2FourPions::t3(Wave4 &q, Wave4 &q1, Wave4 &q2,
                           Wave4 &q3, Wave4 &q4) {
  Wave4 omeQ(q2 + q3 + q4);
  Wave4 rhoQ(q3 + q4);
  double omeS = m2(omeQ);
  double rhoS = m2(rhoQ);

  // Needed to match Herwig++.
  double gM = sqrtpos(rhoM*rhoM - 4*picM*picM) * (rhoM*rhoM - 4*picM*picM)
    / rhoM;
  double dm = (rhoFormFactor1(0) - rhoFormFactor1(rhoM*rhoM) +
                 rhoM*rhoM * rhoFormFactor2(rhoM*rhoM)) / gM;
  return omeW * omeFormFactor(omeS) / (omeD(omeS) * rhoD(rhoS)) *
    pow2(omeM) * (rhoM*rhoM + rhoM*rhoG*dm) *
    ((m2(q,q3) * m2(q1,q4) - m2(q,q4) * m2(q1,q3)) * q2 +
     (m2(q,q4) * m2(q1,q2) - m2(q,q2) * m2(q1,q4)) * q3 +
     (m2(q,q2) * m2(q1,q3) - m2(q,q3) * m2(q1,q2)) * q4);

}

//--------------------------------------------------------------------------

// Return the D function for the a1(1260).

complex HMETau2FourPions::a1D(double s) {

  // rG is defined as the running width.
  double rG = 0;

  // The rho and pion cut off thresholds defined in the fit.
  double piM = 0.16960;
  double rM = 0.83425;

  // Fit of width below three pion mass threshold.
  if (s < piM)
    rG = 0;

  // Fit of width below pion and rho mass threshold.
  else if (s < rM)
    rG = 0.003052*pow3(s - piM)*(1.0 + 151.088*(s - piM) +
                                 174.495*pow2(s - piM));

  // Fit of width above pion and rho mass threshold.
  else
    rG = 2.60817 - 2.47790*s + 0.66539*pow2(s) - 0.0678183*pow3(s) +
        1.66577*(s-1.23701)/s;
  return s - a1M*a1M + complex(0,1) * sqrtpos(s) * rG;

}

//--------------------------------------------------------------------------

// Return the D function for the rho(770).

complex HMETau2FourPions::rhoD(double s) {

  double gQ = sqrtpos(s - 4*picM*picM) * (s - 4*picM*picM) / sqrtpos(s);
  double gM = sqrtpos(rhoM*rhoM - 4*picM*picM) * (rhoM*rhoM - 4*picM*picM)
    / rhoM;
  double dm = (rhoFormFactor1(s) - rhoFormFactor1(rhoM*rhoM) -
                 (s - rhoM*rhoM) * rhoFormFactor2(rhoM*rhoM)) / gM;

  // Ensure complex part is zero below available channel.
  if (s < 4*picM*picM) gQ = 0;
  return s - rhoM*rhoM - rhoM*rhoG*dm + complex(0,1)*rhoM*rhoG*(gQ/gM);

}

//--------------------------------------------------------------------------

// Return the D function for the sigma(800) (just s-wave running width).

complex HMETau2FourPions::sigD(double s) {

  // Sigma decay to two neutral pions for three neutral pion channel.
  double piM = abs(pID[3]) == 111 ? pinM : picM;
  double gQ = sqrtpos(1.0 - 4*piM*piM/s);
  double gM = sqrtpos(1.0 - 4*piM*piM/(sigM*sigM));
  return s - sigM*sigM + complex(0,1)*sigM*sigG*gQ/gM;

}

//--------------------------------------------------------------------------

// Return the D function for the omega(782).

complex HMETau2FourPions::omeD(double s) {

  double g = 0;
  double q = sqrtpos(s);
  double x = q - omeM;

  // Fit of width given in TAUOLA.
  if (s < 1)
    g = 1 + 17.560*x + 141.110*pow2(x) + 894.884*pow3(x) + 4977.35*pow4(x) +
        7610.66*pow5(x) - 42524.4*pow6(x);
  else
    g = -1333.26 + 4860*q - 6000.81*pow2(q) + 2504.97*pow3(q);
  if (g < 0) g = 0;
  return s - omeM*omeM + complex(0,1)*omeM*omeG*g;

}

//--------------------------------------------------------------------------

// Return the form factor for the a1.

double HMETau2FourPions::a1FormFactor(double s) {

  return pow2((1.0 + a1M*a1M/lambda2) / (1.0 + s/lambda2));

}

//--------------------------------------------------------------------------

// Return the form factor for the rho(770) (equivalent to h(s) in TAUOLA).

double HMETau2FourPions::rhoFormFactor1(double s) {

  double f = 0.;
  if (s > 4. * picM * picM) {
    double thr = sqrtpos(1 - 4. * picM * picM / s);
    f = thr * log((1. + thr) / (1. - thr)) * (s - 4. * picM * picM) / M_PI;
  } else if (s < 0.0000001) f = -8. * picM * picM / M_PI;
  return f;

}

//--------------------------------------------------------------------------

// Return the form factor for the rho(770) (equivalent to h(s) derivative).

double HMETau2FourPions::rhoFormFactor2(double s) {

  double f = sqrtpos(1 - 4*picM*picM/s);
  if (s > 4*picM*picM)
    f = f / (M_PI * s) * (s*f + (2*picM*picM + s)*log((1 + f) / (1 - f)));
  else
    f = 0;
  return f;

}

//--------------------------------------------------------------------------

// Return the form factor for the omega(782).

double HMETau2FourPions::omeFormFactor(double /*s*/) {

  return 1.0;

}

//--------------------------------------------------------------------------

// Return the G-functions given in TAUOLA using a piece-wise fit.

double HMETau2FourPions::G(int i, double s) {

  // Break points for the fits.
  double s0(0), s1(0), s2(0), s3(0), s4(0), s5(0);

  // Parameters for the fits.
  double a1(0), b1(0);
  double a2(0), b2(0), c2(0), d2(0), e2(0);
  double a3(0), b3(0), c3(0), d3(0), e3(0);
  double a4(0), b4(0);
  double a5(0), b5(0);

  // Three neutral pion parameters.
  if (i == 1) {
    s0 = 0.614403;      s1 = 0.656264;  s2 = 1.57896;
    s3 = 3.08198;       s4 = 3.12825;   s5 = 3.17488;
    a1 = -23383.7;      b1 = 38059.2;
    a2 = 230.368;       b2 = -4.39368;  c2 = 687.002;
    d2 = -732.581;      e2 = 207.087;
    a3 = 1633.92;       b3 = -2596.21;  c3 = 1703.08;
    d3 = -501.407;      e3 = 54.5919;
    a4 = -2982.44;      b4 = 986.009;
    a5 = 6948.99;       b5 = -2188.74;
  }

  // Three charged pion parameters.
  else if (i == 2) {
    s0 = 0.614403;      s1 = 0.635161;  s2 = 2.30794;
    s3 = 3.08198;       s4 = 3.12825;   s5 = 3.17488;
    a1 = -54171.5;      b1 = 88169.3;
    a2 = 454.638;       b2 = -3.07152;  c2 = -48.7086;
    d2 = 81.9702;       e2 = -24.0564;
    a3 = -162.421;      b3 = 308.977;   c3 = -27.7887;
    d3 = -48.5957;      e3 = 10.6168;
    a4 = -2650.29;      b4 = 879.776;
    a5 = 6936.99;       b5 = -2184.97;
  }

  // Omega mediated three charged pion parameters.
  else if (i == 3) {
    s0 = 0.81364;       s1 = 0.861709;  s2 = 1.92621;
    s3 = 3.08198;       s4 = 3.12825;   s5 = 3.17488;
    a1 = -84888.9;      b1 = 104332;
    a2 = 2698.15;       b2 = -3.08302;  c2 = 1936.11;
    d2 = -1254.59;      e2 = 201.291;
    a3 = 7171.65;       b3 = -6387.9;   c3 = 3056.27;
    d3 = -888.63;       e3 = 108.632;
    a4 = -5607.48;      b4 = 1917.27;
    a5 = 26573; b5 = -8369.76;
  }

  // Return the appropriate fit.
  if (s < s0)
    return 0.0;
  else if (s < s1)
   return a1 + b1*s;
  else if (s < s2)
    return a2*pow(s,b2) + c2*pow2(s) + d2*pow3(s) + e2*pow4(s);
  else if (s < s3)
    return a3 + b3*s + c3*pow2(s) + d3*pow3(s) + e3*pow4(s);
  else if (s < s4)
    return a4 + b4*s;
  else if (s < s5)
    return a5 + b5*s;
  else
    return 0.0;

}

//==========================================================================

// Tau decay matrix element for tau decay into five pions using the model given
// in hep-ph/0602162v1.

// Ja(q,q1,q2,q3,q4,q5): current through rho and omega resonances
// Jb(q,q1,q2,q3,q4,q5): current through a1 and sigma resonances
// breitWigner(s,M,G): s-wave Breit-Wigner assuming massless products
// a1M: on-shell mass of the a1 resonance
// a1G: width of the a1 resonance
// rhoX: mass and width of the rho resonance
// omegaX: mass, width, and weight of the omega resonance
// sigmaX: mass, width, and weight of the sigma resonance

//--------------------------------------------------------------------------

// Initialize constants for the helicity matrix element.

void HMETau2FivePions::initConstants() {

  // pi-, pi-, pi+, pi+, pi- decay.
  if (abs(pID[2]) == 211 && abs(pID[3]) == 211 && abs(pID[4]) == 211 &&
      abs(pID[5]) == 211 && abs(pID[6]) == 211)
    DECAYWEIGHTMAX = 4e4;
  // pi+, pi-, pi0, pi-, pi0 decay.
  else if (abs(pID[2]) == 111 && abs(pID[3]) == 111 && abs(pID[4]) == 211 &&
           abs(pID[5]) == 211 && abs(pID[6]) == 211)
    DECAYWEIGHTMAX = 1e7;
  // pi0, pi0, pi-, pi0, pi0 decay.
  else if (abs(pID[2]) == 111 && abs(pID[3]) == 111 && abs(pID[4]) == 111 &&
           abs(pID[5]) == 111 && abs(pID[6]) == 211)
    DECAYWEIGHTMAX = 1e5;

  // Set resonances.
  a1M    = 1.260; a1G    = 0.400;
  rhoM   = 0.776; rhoG   = 0.150;
  omegaM = 0.782; omegaG = 0.0085; omegaW = 11.5;
  sigmaM = 0.800; sigmaG = 0.600;  sigmaW = 1;

}

//--------------------------------------------------------------------------

// Initialize the hadronic current for the helicity matrix element.

void HMETau2FivePions::initHadronicCurrent(vector<HelicityParticle>& p) {

  vector< Wave4 > u2;

  Wave4 q(p[2].p() + p[3].p() + p[4].p() + p[5].p() + p[6].p());
  Wave4 q2(p[2].p()), q3(p[3].p()), q4(p[4].p()), q5(p[5].p()), q6(p[6].p());
  // pi-, pi-, pi+, pi+, pi- decay.
  if (abs(pID[2]) == 211 && abs(pID[3]) == 211 && abs(pID[4]) == 211 &&
      abs(pID[5]) == 211 && abs(pID[6]) == 211)
    u2.push_back(Jb(q, q2, q3, q5, q6, q4) + Jb(q, q4, q3, q5, q6, q2)
                 + Jb(q, q2, q4, q5, q6, q3) + Jb(q, q2, q3, q6, q5, q4)
                 + Jb(q, q4, q3, q6, q5, q2) + Jb(q, q2, q4, q6, q5, q3));
  // pi+, pi-, pi0, pi-, pi0 decay.
  else if (abs(pID[2]) == 111 && abs(pID[3]) == 111 && abs(pID[4]) == 211 &&
           abs(pID[5]) == 211 && abs(pID[6]) == 211)
    u2.push_back(Ja(q, q6, q4, q2, q5, q3) + Ja(q, q6, q5, q2, q4, q3)
                 + Ja(q, q6, q4, q3, q5, q2) + Ja(q, q6, q5, q3, q4, q2)
                 + Jb(q, q4, q5, q6, q2, q3) + Jb(q, q2, q3, q4, q6, q5)
                 + Jb(q, q2, q3, q5, q6, q4));
  // pi0, pi0, pi-, pi0, pi0 decay.
  else if (abs(pID[2]) == 111 && abs(pID[3]) == 111 && abs(pID[4]) == 111 &&
           abs(pID[5]) == 111 && abs(pID[6]) == 211)
    u2.push_back(Jb(q, q2, q3, q6, q4, q5) + Jb(q, q5, q3, q6, q4, q2)
                 + Jb(q, q3, q4, q6, q2, q5) + Jb(q, q2, q4, q6, q3, q5)
                 + Jb(q, q2, q5, q6, q4, q3) + Jb(q, q4, q5, q6, q2, q3));

  u.push_back(u2);

}

//--------------------------------------------------------------------------

// Return the omega-rho hadronic current.

Wave4 HMETau2FivePions::Ja(Wave4 &q, Wave4 &q1, Wave4 &q2,
                           Wave4 &q3, Wave4 &q4, Wave4 &q5) {

  Wave4 j = epsilon(q1, q2, q3);
  return omegaW * (breitWigner(m2(q), a1M, a1G)
                   * breitWigner(m2(q1 + q2 + q3), omegaM, omegaG)
                   * breitWigner(m2(q4 + q5), rhoM, rhoG)
                   * epsilon(q4 - q5, j, q)
                   * (breitWigner(m2(q2 + q3), rhoM, rhoG)
                      + breitWigner(m2(q1 + q3), rhoM, rhoG)
                      + breitWigner(m2(q1 + q2), rhoM, rhoG)));

}

//--------------------------------------------------------------------------

// Return the a1-sigma hadronic current.

Wave4 HMETau2FivePions::Jb(Wave4 &q, Wave4 &q1, Wave4 &q2,
                           Wave4 &q3, Wave4 &q4, Wave4 &q5) {

  double s = m2(q);
  Wave4  a1Q = q1 + q2 + q3;
  double a1S = m2(a1Q);
  Wave4 j = (m2(q2, q1 - q3) / a1S * a1Q - q1 + q3)
    * breitWigner(m2(q1 + q3), rhoM, rhoG)
    + (m2(q1, q2 - q3) / a1S * a1Q - q2 + q3)
    * breitWigner(m2(q2 + q3), rhoM, rhoG);
  j = (j * gamma[4] * q / s) * q - j;
  return sigmaW * (breitWigner(s, a1M, a1G) * breitWigner(a1S, a1M, a1G)
                   * breitWigner(m2(q4 + q5), sigmaM, sigmaG) * j);

}

  complex HMETau2FivePions::breitWigner(double s, double M, double G) {

    return M * M / (M * M - s - complex(0, 1) * M * G);

}

//==========================================================================

} // end namespace Pythia8
