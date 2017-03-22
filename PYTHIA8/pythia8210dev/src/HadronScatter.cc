// HadronScatter.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia8/HadronScatter.h"

namespace Pythia8 {

//==========================================================================

// The SigmaPartialWave class
//  Reads in tables of partial wave data to provide dSigma/dCos(theta)
//  The generic classes of process are:
//    process = 0 (pi-pi), 1 (pi-K), 2 (pi-N)
//  Subprocesses are defined (along with isospin coefficients) in:
//    setupSubprocesses();
//  Individual subprocesses are selected using:
//    setSubprocess(subprocess); or setSubprocess(PDG1, PDG2);
//  Internally, there are two std::map's, to convert between:
//    subprocess <==> PDG1, PDG2
//
//  Data are read in from files:
//   Lines starting with a '#' are comments
//   Lines starting with 'set' provide options:
//    set eType   [Wcm | Tlab] - energy bins in Wcm or Tlab
//    set eUnit   [MeV | GeV]  - energy unit
//    set input   [eta,delta | Sn,delta  | Tr,Ti | mod,phi ]
//                             - format of columns in partial waves
//    set dUnit   [deg | rad]  - unit of phase shifts
//    set norm    [0 | 1]      - normalisation
//   Column headers give L,2I[,2J] (2J for e.g. piN)
//   Input types: Sn,delta -> Sn = 1 - eta^2
//                mod,phi  -> amplitude T_L = |T_L| exp(i phi_L)
//   Normalisation: 0 -> dSigma/dOmega = 1 / k^2 |T_L|^2
//                  1 -> dSigma/dOmega = 16 / s  |T_L|^2
//
//  Internally data is stored as (J = 0 for spinless):
//   pwData[L * LSHIFT + 2I * ISHIFT + J][energy_bin_centre] = T
//  where the energy is Wcm in GeV.
//
//  This is stored using std::map's, to take into account that not all
//  L,I,J states are always present (e.g. negligable contributions or
//  conservation rules) and that bin sizes are not fixed.
//
//  Re[T] and Im[T] are interpolated between bins and extrapolated down to
//  threshold from the first two bins. Above energy_bin_centre of the final
//  bin, no extrapolation is done and the final bin value is always used.
//
//  A simple scheme to provide correct distributions for cos(theta) at a
//  given CM energy is included. Efficiency is not too bad, but can likely
//  be greatly improved.
//
//  For each subprocess, a grid in bins of Wcm and cos(theta) is setup with:
//    setupGrid();
//  The size of the grid is set by the constants:
//    const double SigmaPartialWave::WCMBIN;
//    const double SigmaPartialWave::CTBIN;
//  For each bin of (Wcm, ct), the maximum sigma elastic is found by
//  splitting this bin into subbins multiple times, controlled by:
//    const int SigmaPartialWave::SUBBIN;
//    const int SigmaPartialWave::ITER
//  With the final maxium sigma elastic given by this value multipled by
//  a safety factor:
//    const double SigmaPartialWave::GRIDSAFETY
//
//  To pick values of cos(theta) for a given CM energy, a:
//    pickCosTheta(Wcm);
//  function is provided. The above grid is used as an overestimate, to
//  pick properly distributed values of cos(theta).

//--------------------------------------------------------------------------

// Constants

// pwData[L * LSHIFT + 2I * ISHIFT + J]
const int     SigmaPartialWave::LSHIFT     = 1000000;
const int     SigmaPartialWave::ISHIFT     = 1000;

// Convert GeV^-2 to mb
const double  SigmaPartialWave::CONVERT2MB = 0.389380;

// Size of bin in Wcm and cos(theta)
const double  SigmaPartialWave::WCMBIN     = 0.005;
const double  SigmaPartialWave::CTBIN      = 0.2;
// Number of subbins and iterations
const int     SigmaPartialWave::SUBBIN     = 2;
const int     SigmaPartialWave::ITER       = 2;
// Safety value to add on to grid maxima
const double  SigmaPartialWave::MASSSAFETY = 0.001;
const double  SigmaPartialWave::GRIDSAFETY = 0.05;


//--------------------------------------------------------------------------

// Perform initialization and store pointers.

bool SigmaPartialWave::init(int processIn, string xmlPath, string filename,
                            Info *infoPtrIn, ParticleData *particleDataPtrIn,
                            Rndm *rndmPtrIn) {
  // Store incoming pointers
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;

  // Check incoming process is okay
  if (processIn < 0 || processIn > 2) {
    infoPtr->errorMsg("Error in SigmaPartialWave::init: "
      "unknown process");
    return false;
  }
  process = processIn;

  // Setup subprocesses and isospin coefficients
  setupSubprocesses();
  setSubprocess(0);

  // Read in partial-wave data
  if (!readFile(xmlPath, filename)) return false;

  // Setup vector for Legendre polynomials
  PlVec.resize(Lmax);
  if (Lmax > 0) PlVec[0] = 1.;
  // And derivatives if needed
  if (process == 2) {
    PlpVec.resize(Lmax);
    if (Lmax > 0) PlpVec[0] = 0.;
    if (Lmax > 1) PlpVec[1] = 1.;
  }

  // Setup grid for integration
  setupGrid();

  return true;
}


//--------------------------------------------------------------------------

// Read input data file

bool SigmaPartialWave::readFile(string xmlPath, string filename) {
  // Create full path and open file
  string fullPath = xmlPath + filename;
  ifstream ifs(fullPath.c_str());
  if (!ifs.good()) {
    infoPtr->errorMsg("Error in SigmaPartialWave::init: "
      "could not read data file");
    return false;
  }

  // Default unit settings
  int eType = 0;  // 0 = Wcm, 1 = Tlab
  int eUnit = 0;  // 0 = GeV, 1 = MeV
  int input = 0;  // 0 = eta, delta; 1 = Sn, delta (Sn = 1 - eta^2);
                  // 2 = Treal, Tim, 3 = mod, phi
  int dUnit = 0;  // 0 = deg, 1 = rad
      norm  = 0;  // 0 = standard, 1 = sqrt(s) / sqrt(s - 4Mpi^2)

  // Once we have a header line, each column corresponds to
  // values of L, I and J.
  Lmax = Imax = 0;
  binMax = 0.;
  vector < int > Lvec, Ivec, Jvec;

  // Parse the file
  string line;
  while (ifs.good()) {
    // Get line, convert to lowercase and strip leading whitespace
    getline(ifs, line);
    for (unsigned int i = 0; i < line.length(); i++)
      line[i] = tolower(line[i]);
    string::size_type startPos = line.find_first_not_of("  ");
    if (startPos != string::npos) line = line.substr(startPos);
    // Skip blank lines and lines that start with '#'
    if (line.length() == 0 || line[0] == '#') continue;

    // Tokenise line on whitespace (spaces or tabs)
    string lineT = line;
    vector < string > token;
    while (true) {
      startPos = lineT.find_first_of("   ");
      token.push_back(lineT.substr(0, startPos));
      if (startPos == string::npos) break;
      startPos = lineT.find_first_not_of("       ", startPos + 1);
      if (startPos == string::npos) break;
      lineT = lineT.substr(startPos);
    }

    // Settings
    if (token[0] == "set") {
      bool badSetting = false;

      // eType
      if        (token[1] == "etype") {
        if      (token[2] == "wcm")       eType = 0;
        else if (token[2] == "tlab")      eType = 1;
        else    badSetting = true;

      // eUnit
      } else if (token[1] == "eunit") {
        if      (token[2] == "gev")       eUnit = 0;
        else if (token[2] == "mev")       eUnit = 1;
        else    badSetting = true;

      // input
      } else if (token[1] == "input") {
        if      (token[2] == "eta,delta") input = 0;
        else if (token[2] == "sn,delta")  input = 1;
        else if (token[2] == "tr,ti")     input = 2;
        else if (token[2] == "mod,phi")   input = 3;
        else    badSetting = true;

      // dUnit
      } else if (token[1] == "dunit") {
        if      (token[2] == "deg")       dUnit = 0;
        else if (token[2] == "rad")       dUnit = 1;
        else    badSetting = true;

      // norm
      } else if (token[1] == "norm") {
        if      (token[2] == "0")         norm = 0;
        else if (token[2] == "1")         norm = 1;
        else    badSetting = true;
      }

      // Bad setting
      if (badSetting) {
        infoPtr->errorMsg("Error in SigmaPartialWave::init: "
          "bad setting line");
        return false;
      }
      continue;
    }

    // Header line
    if (line.substr(0, 1).find_first_of("0123456789.") != 0) {
      // Clear current stored L,2I,2J values
      Lvec.clear(); Ivec.clear(); Jvec.clear();

      // Parse header
      bool badHeader = false;
      for (unsigned int i = 1; i < token.size(); i++) {
        // Extract L
        startPos = token[i].find_first_of(",");
        if (startPos == string::npos) { badHeader = true; break; }
        string Lstr = token[i].substr(0, startPos);
        token[i] = token[i].substr(startPos + 1);
        // Extract 2I
        string Istr;
        startPos = token[i].find_first_of(",     ");
        if (startPos == string::npos) {
          Istr = token[i];
          token[i] = "";
        } else {
          Istr = token[i].substr(0, startPos);
          token[i] = token[i].substr(startPos + 1);
        }
        // Extract 2J
        string Jstr("0");
        if (token[i].length() != 0) Jstr = token[i];

        // Convert to integers and store
        int L, I, J;
        stringstream Lss(Lstr); Lss >> L;
        stringstream Iss(Istr); Iss >> I;
        stringstream Jss(Jstr); Jss >> J;
        if (Lss.fail() || Iss.fail() || Jss.fail())
          { badHeader = true; break; }
        Lvec.push_back(L);
        Ivec.push_back(I);
        Jvec.push_back(J);
        Lmax = max(Lmax, L);
        Imax = max(Imax, I);
      }
      if (badHeader) {
        infoPtr->errorMsg("Error in SigmaPartialWave::init: "
          "malformed header line");
        return false;
      }

    // Data line
    } else {
      bool badData = false;

      // Check there are the correct number of columns
      if (token.size() != 2 * Lvec.size() + 1) badData = true;

      // Extract energy
      double eNow = 0.;
      if (!badData) {
        stringstream eSS(token[0]);
        eSS >> eNow;
        if (eSS.fail()) badData = true;
        // Convert to GeV if needed
        if (eUnit == 1) eNow *= 1e-3;
        // Convert to Wcm if needed
        if (eType == 1) eNow = sqrt(2. * mB * eNow + pow2(mA + mB));
        binMax = max(binMax, eNow);
      }

      // Extract eta/phase shifts
      if (!badData) {
        for (unsigned int i = 1; i < token.size(); i += 2) {
          // L,2I,2J
          int LIJidx = (i - 1) / 2;
          int L = Lvec[LIJidx];
          int I = Ivec[LIJidx];
          int J = Jvec[LIJidx];

          double i1, i2;
          stringstream i1SS(token[i]);     i1SS >> i1;
          stringstream i2SS(token[i + 1]); i2SS >> i2;
          if (i1SS.fail() || i2SS.fail()) { badData = true; break; }

          // Sn to eta
          if (input == 1) i1 = sqrt(1. - i1);
          // Degrees to radians
          if ((input == 0 || input == 1 || input == 3) &&
              dUnit == 0) i2 *= M_PI / 180.;

          // Convert to Treal and Timg
          complex T(0., 0.);
          if (input == 0 || input == 1) {
            T = (i1 * exp(2. * complex(0., 1.) * i2) - 1.) /
                2. / complex(0., 1.);
          } else if (input == 2) {
            T = complex(i1, i2);
          } else if (input == 3) {
            T = i1 * exp(complex(0., 1.) * i2);
          }

          // Store
          pwData[L * LSHIFT + I * ISHIFT + J][eNow] = T;
        }
      }
      if (badData) {
        infoPtr->errorMsg("Error in SigmaPartialWave::init: "
          "malformed data line");
        return false;
      }
    }
  }

  // Make sure it was EOF that caused us to end
  if (!ifs.eof()) { ifs.close(); return false; }

  // Maximum values of L and I
  Lmax++; Imax++;

  return true;
}


//--------------------------------------------------------------------------

// Setup isospin coefficients and subprocess mapping

void SigmaPartialWave::setupSubprocesses() {

  // Setup isospin coefficients
  switch (process) {
  // pi-pi
  case 0:
    // Map subprocess to incoming
    subprocessMax = 6;
    sp2in[0] = pair < int, int > ( 211,  211);
    sp2in[1] = pair < int, int > ( 211, -211);
    sp2in[2] = pair < int, int > ( 211,  111);
    sp2in[3] = pair < int, int > ( 111,  111);
    sp2in[4] = pair < int, int > (-211,  111);
    sp2in[5] = pair < int, int > (-211, -211);
    // Incoming to subprocess
    for (int i = 0; i < subprocessMax; i++) in2sp[sp2in[i]] = i;
    // Isospin coefficients
    isoCoeff[0][0] = 0.;    isoCoeff[0][2] = 0.;    isoCoeff[0][4] = 1.;
    isoCoeff[1][0] = 1./3.; isoCoeff[1][2] = 1./2.; isoCoeff[1][4] = 1./6.;
    isoCoeff[2][0] = 0.;    isoCoeff[2][2] = 1./2.; isoCoeff[2][4] = 1./2.;
    isoCoeff[3][0] = 1./3.; isoCoeff[3][2] = 0.;    isoCoeff[3][4] = 2./3.;
    isoCoeff[4][0] = 0.;    isoCoeff[4][2] = 1./2.; isoCoeff[4][4] = 1./2.;
    isoCoeff[5][0] = 0.;    isoCoeff[5][2] = 0.;    isoCoeff[5][4] = 1.;

    break;

  // pi-K and pi-N
  case 1: case 2:
    int id1, id2;
    if (process == 1) { id1 = 321;  id2 = 311;  }
    else              { id1 = 2212; id2 = 2112; }

    // Map subprocess to incoming
    subprocessMax = 12;
    sp2in[0] = pair < int, int > ( 211, id1);
    sp2in[1] = pair < int, int > ( 211, id2);
    sp2in[2] = pair < int, int > ( 111, id1);
    sp2in[3] = pair < int, int > ( 111, id2);
    sp2in[4] = pair < int, int > (-211, id1);
    sp2in[5] = pair < int, int > (-211, id2);
    // Isospin coefficients
    isoCoeff[0][1]  = 0.;      isoCoeff[0][3]  = 1.;
    isoCoeff[1][1]  = 2. / 3.; isoCoeff[1][3]  = 1. / 3.;
    isoCoeff[2][1]  = 1. / 3.; isoCoeff[2][3]  = 2. / 3.;
    isoCoeff[3][1]  = 1. / 3.; isoCoeff[3][3]  = 2. / 3.;
    isoCoeff[4][1]  = 2. / 3.; isoCoeff[4][3]  = 1. / 3.;
    isoCoeff[5][1]  = 0.;      isoCoeff[5][3]  = 1.;
    // Antiparticles
    for (int i = 0; i < 6; i++) {
      id1 = ((sp2in[i].first == 111) ? +1 : -1) * sp2in[i].first;
      sp2in[i + 6] = pair < int, int > (id1, -sp2in[i].second);
      isoCoeff[i + 6] = isoCoeff[i];
    }
    // Map incoming to subprocess
    for (int i = 0; i < subprocessMax; i++) in2sp[sp2in[i]] = i;

    break;
  }

  return;
}


//--------------------------------------------------------------------------

// Setup grids for integration

void SigmaPartialWave::setupGrid() {
  // Reset sigma maximum
  sigElMax = 0.;

  // Go through each subprocess
  gridMax.resize(subprocessMax);
  gridNorm.resize(subprocessMax);
  for (int sp = 0; sp < subprocessMax; sp++) {
    // Setup subprocess
    setSubprocess(sp);

    // Bins in Wcm
    int nBin1 = int( (binMax - mA - mB) / WCMBIN );
    gridMax[subprocess].resize(nBin1);
    gridNorm[subprocess].resize(nBin1);
    for (int n1 = 0; n1 < nBin1; n1++) {
      // Bin lower and upper
      double bl1 = mA + mB + double(n1) * WCMBIN;
      double bu1 = bl1 + WCMBIN;

      // Bins in cos(theta)
      int    nBin2 = int( 2. / CTBIN );
      gridMax[subprocess][n1].resize(nBin2);
      for (int n2 = 0; n2 < nBin2; n2++) {
        // Bin lower and upper
        double bl2 = -1. + double(n2) * CTBIN;
        double bu2 = bl2 + CTBIN;

        // Find maximum
        double maxSig = 0.;
        double bl3 = bl1, bu3 = bu1, bl4 = bl2, bu4 = bu2;
        for (int iter = 0; iter < ITER; iter++) {
          int    i3Save = -1, i4Save = -1;
          double step3 = (bu3 - bl3) / double(SUBBIN);
          double step4 = (bu4 - bl4) / double(SUBBIN);
          for (int i3 = 0; i3 <= SUBBIN; i3++) {
            double Wcm = bl3 + double(i3) * step3;
            for (int i4 = 0; i4 <= SUBBIN; i4++) {
              double ct = bl4 + double(i4) * step4;
              double ds = dSigma(Wcm, ct);
              if (ds > maxSig) {
                i3Save = i3;
                i4Save = i4;
                maxSig = ds;
              }
            }
          }
          // Set new min/max
          if (i3Save == -1 && i4Save == -1) break;
          if (i3Save > -1) {
            bl3 = bl3 + ((i3Save == 0)      ? 0. : i3Save - 1.) * step3;
            bu3 = bl3 + ((i3Save == SUBBIN) ? 1. : 2.)          * step3;
          }
          if (i4Save > -1) {
            bl4 = bl4 + ((i4Save == 0)      ? 0. : i4Save - 1.) * step4;
            bu4 = bl4 + ((i4Save == SUBBIN) ? 1. : 2.)          * step4;
          }
        } // for (iter)

        // Save maximum value
        gridMax[subprocess][n1][n2]  = maxSig * (1. + GRIDSAFETY);
        gridNorm[subprocess][n1]    += maxSig * (1. + GRIDSAFETY) * CTBIN;
        sigElMax = max(sigElMax, maxSig);

      } // for (n2)
    } // for (n1)
  } // for (sp)

  return;
}


//--------------------------------------------------------------------------

// Pick a cos(theta) value

double SigmaPartialWave::pickCosTheta(double Wcm) {
  // Find grid bin in Wcm
  int WcmBin = int((Wcm - mA - mB) / WCMBIN);
  if (WcmBin < 0) WcmBin = 0;
  if (WcmBin >= int(gridMax[subprocess].size()))
    WcmBin = int(gridMax[subprocess].size() - 1);

  // Pick a value of cos(theta)
  double ct, wgt;

  do {
    // Sample from overestimate and inverse
    double y   = rndmPtr->flat() * gridNorm[subprocess][WcmBin];
    double sum = 0.;
    int    ctBin;
    for (ctBin = 0; ctBin < int(2. / CTBIN); ctBin++) {
      if (sum + CTBIN * gridMax[subprocess][WcmBin][ctBin] > y) break;
      sum += CTBIN * gridMax[subprocess][WcmBin][ctBin];
    }

    // Linear interpolation
    double x1 = -1. + CTBIN * double(ctBin);
    double y1 = sum;
    double x2 = x1 + CTBIN;
    double y2 = sum + CTBIN * gridMax[subprocess][WcmBin][ctBin];
           ct = (x2 - x1) / (y2 - y1) * (y - y1) + x1;
    wgt = dSigma(Wcm, ct) / gridMax[subprocess][WcmBin][ctBin];
    if (wgt >= 1.) {
      infoPtr->errorMsg("Warning in SigmaPartialWave::pickCosTheta: "
        "weight above unity");
      break;
    }
  } while (wgt <= rndmPtr->flat());

  return ct;
}

//--------------------------------------------------------------------------

// Set subprocess

bool SigmaPartialWave::setSubprocess(int spIn) {
  if (sp2in.find(spIn) == sp2in.end()) return false;
  subprocess = spIn;
  pair < int, int > in = sp2in[spIn];
  idA = in.first;
  mA  = particleDataPtr->m0(idA);
  idB = in.second;
  mB  = particleDataPtr->m0(idB);
  return true;
}

bool SigmaPartialWave::setSubprocess(int idAin, int idBin) {
  pair < int, int > in = pair < int, int > (idAin, idBin);
  if (in2sp.find(in) == in2sp.end()) {
    // Try the other way around as well
    swap(in.first, in.second);
    if (in2sp.find(in) == in2sp.end()) return false;
  }
  subprocess = in2sp[in];
  idA = idAin;
  mA  = particleDataPtr->m0(idA);
  idB = idBin;
  mB  = particleDataPtr->m0(idB);
  return true;
}

//--------------------------------------------------------------------------

// Calculate: mode = 0 (sigma elastic), 1 (sigma total), 2 (dSigma/dcTheta)

double SigmaPartialWave::sigma(int mode, double Wcm, double cTheta) {
  // Below threshold, return 0
  if (Wcm < (mA + mB + MASSSAFETY)) return 0.;

  // Return values
  complex amp[2] = { complex(0., 0.) };
  double  sig    = 0.;

  // Kinematic variables
  double s  = pow2(Wcm);
  double k2 = (s - pow2(mB + mA)) * (s - pow2(mB - mA)) / 4. / s;

  // Precompute all required Pl and Pl' values
  double sTheta = 0.;
  if (mode == 2) {
    if (process == 2) sTheta = sqrt(1. - pow2(cTheta));
    legendreP(cTheta, ((process == 2) ? true : false));
  }

  // Loop over L
  for (int L = 0; L < Lmax; L++) {

    // Loop over J (only J = 0 for spinless)
    complex ampJ[2] = { complex(0., 0.) };
    int Jstart = (process != 2) ? 0 : 2 * L - 1;
    int Jend   = (process != 2) ? 1 : 2 * L + 2;
    int Jstep  = (process != 2) ? 1 : 2;
    int Jcount = 0;
    for (int J = Jstart; J < Jend; J += Jstep, Jcount++) {

      // Loop over isospin coefficients
      for (int I = 0; I < Imax; I++) {
        if (isoCoeff[subprocess][I] == 0.) continue;

        // Check wave exists
        int LIJ = L * LSHIFT + I * ISHIFT + J;
        if (pwData.find(LIJ) == pwData.end()) continue;

        // Extrapolation / interpolation (not for last bin)
        map < double, complex >::iterator it = pwData[LIJ].upper_bound(Wcm);
        if (it == pwData[LIJ].begin()) ++it;
        double ar, ai;
        if (it == pwData[LIJ].end()) {
          ar = (--it)->second.real();
          ai = it->second.imag();
        } else {
          double  eA   = it->first;
          complex ampA = (it--)->second;
          double  eB   = it->first;
          complex ampB = it->second;

          ar = (ampA.real() - ampB.real()) / (eA - eB) *
               (Wcm - eB) + ampB.real();
          ai = (ampA.imag() - ampB.imag()) / (eA - eB) *
               (Wcm - eB) + ampB.imag();
        }

        // Isospin sum
        ampJ[Jcount] += isoCoeff[subprocess][I] * complex(ar, ai);
      }
    }

    // Partial wave sum. Sigma elastic
    if (mode == 0) {
      if        (process == 0 || process == 1) {
        sig += (2. * L + 1.) * (ampJ[0] * conj(ampJ[0])).real();
      } else if (process == 2) {
        sig += ( (L + 0.) * (ampJ[0] * conj(ampJ[0])).real() +
                 (L + 1.) * (ampJ[1] * conj(ampJ[1])).real() );
      }

    // Sigma total
    } else if (mode == 1) {
      if        (process == 0 || process == 1) {
        sig += (2. * L + 1.) * ampJ[0].imag();
      } else if (process == 2) {
        sig += ( (L + 0.) * ampJ[0].imag() + (L + 1.) * ampJ[1].imag() );
      }

    // dSigma
    } else if (mode == 2) {
      if        (process == 0 || process == 1) {
        amp[0] += (2. * L + 1.) * ampJ[0] * PlVec[L];
      } else if (process == 2) {
        amp[0] += ((L + 0.) * ampJ[0] + double(L + 1.) * ampJ[1]) * PlVec[L];
        amp[1] += complex(0., 1.) * (ampJ[1] - ampJ[0]) * sTheta * PlpVec[L];
      }
    }

  } // for (L)

  // Normalisation and return
  if (mode == 0 || mode == 1) {
    if      (norm == 0)  sig *= 4.  * M_PI / k2 * CONVERT2MB;
    else if (norm == 1)  sig *= 64. * M_PI / s  * CONVERT2MB;

  } else if (mode == 2) {
    sig = (amp[0] * conj(amp[0])).real() + (amp[1] * conj(amp[1])).real();
    if      (norm == 0) sig *= 2.  * M_PI / k2 * CONVERT2MB;
    else if (norm == 1) sig *= 32. * M_PI / s  * CONVERT2MB;
  }
  // Half for identical
  return ((idA == idB) ? 0.5 : 1.) * sig;
}


//--------------------------------------------------------------------------

// Bonnet's recursion formula for Legendre polynomials and derivatives

void SigmaPartialWave::legendreP(double ct, bool deriv) {
  if (Lmax > 1) PlVec[1] = ct;
  for (int L = 2; L < Lmax; L++) {
    PlVec[L] = ( (2. * L - 1.) * ct * PlVec[L - 1] -
                 (L - 1.) * PlVec[L - 2] ) / double(L);
    if (deriv)
      PlpVec[L] = ( (2. * L - 1.) * (PlVec[L - 1] + ct * PlpVec[L - 1]) -
                    (L - 1.) * PlpVec[L - 2] ) / double(L);
  }
  return;
}


//==========================================================================

// HadronScatter class

//--------------------------------------------------------------------------

// Perform initialization and store pointers.

bool HadronScatter::init(Info* infoPtrIn, Settings& settings,
                         Rndm *rndmPtrIn, ParticleData *particleDataPtr) {
  // Save incoming pointers
  infoPtr = infoPtrIn;
  rndmPtr = rndmPtrIn;

  // Main settings
  doHadronScatter = settings.flag("HadronScatter:scatter");
  afterDecay      = settings.flag("HadronScatter:afterDecay");
  allowDecayProd  = settings.flag("HadronScatter:allowDecayProd");
  scatterRepeat   = settings.flag("HadronScatter:scatterRepeat");
  // Hadron selection
  hadronSelect    = settings.mode("HadronScatter:hadronSelect");
  Npar            = settings.parm("HadronScatter:N");
  kPar            = settings.parm("HadronScatter:k");
  pPar            = settings.parm("HadronScatter:p");
  // Scattering probability
  scatterProb     = settings.mode("HadronScatter:scatterProb");
  jPar            = settings.parm("HadronScatter:j");
  rMax            = settings.parm("HadronScatter:rMax");
  rMax2           = rMax * rMax;
  doTile          = settings.flag("HadronScatter:tile");

  // String fragmentation and MPI settings
  pTsigma         = 2.0 * settings.parm("StringPT:sigma");
  pTsigma2        = pTsigma * pTsigma;
  double pT0ref   = settings.parm("MultipartonInteractions:pT0ref");
  double eCMref   = settings.parm("MultipartonInteractions:eCMref");
  double eCMpow   = settings.parm("MultipartonInteractions:eCMpow");
  double eCMnow   = infoPtr->eCM();
  pT0MPI          = pT0ref * pow(eCMnow / eCMref, eCMpow);

  // Tiling
  double mp2 = particleDataPtr->m0(111) * particleDataPtr->m0(111);
  double eA  = infoPtr->eA();
  double eB  = infoPtr->eB();
  double pzA =  sqrt(eA * eA - mp2);
  double pzB = -sqrt(eB * eB - mp2);
  yMax = 0.5 * log((eA + pzA) / (eA - pzA));
  yMin = 0.5 * log((eB + pzB) / (eB - pzB));
  // Size in y and phi
  if (doTile) {
    ytMax  = int((yMax - yMin) / rMax);
    ytSize = (yMax - yMin) / double(ytMax);
    ptMax  = int(2. * M_PI / rMax);
    ptSize = 2. * M_PI / double(ptMax);
  } else {
    ytMax  = 1;
    ytSize = yMax - yMin;
    ptMax  = 1;
    ptSize = 2. * M_PI;
  }
  // Initialise tiles
  tile.resize(ytMax);
  for (int yt = 0; yt < ytMax; yt++) tile[yt].resize(ptMax);

  // Find path to data files, i.e. xmldoc directory location.
  // Environment variable takes precedence, else use constructor input.
  // XXX - as in Pythia.cc, but not passed around in e.g. Info/Settings,
  //       so redo here
  string xmlPath = "";
  const char* PYTHIA8DATA = "PYTHIA8DATA";
  char* envPath = getenv(PYTHIA8DATA);
  if (envPath != 0 && *envPath != '\0') {
    int i = 0;
    while (*(envPath+i) != '\0') xmlPath += *(envPath+(i++));
  } else xmlPath = "../xmldoc";
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";

  // Hadron scattering partial wave cross sections
  if ( !sigmaPW[0].init(0, xmlPath, "pipi-Froggatt.dat",
                        infoPtr, particleDataPtr, rndmPtr) ) return false;
  if ( !sigmaPW[1].init(1, xmlPath, "piK-Estabrooks.dat",
                        infoPtr, particleDataPtr, rndmPtr) ) return false;
  if ( !sigmaPW[2].init(2, xmlPath, "piN-SAID-WI08.dat",
                        infoPtr, particleDataPtr, rndmPtr) ) return false;
  sigElMax = 0.;
  sigElMax = max(sigElMax, sigmaPW[0].getSigmaElMax());
  sigElMax = max(sigElMax, sigmaPW[1].getSigmaElMax());
  sigElMax = max(sigElMax, sigmaPW[2].getSigmaElMax());

  // DEBUG
  debugOutput();

  return true;
}


//--------------------------------------------------------------------------

// Debug output

void HadronScatter::debugOutput() {
  // Print settings
  cout << "Hadron scattering:" << endl
       << " scatter        = " << ((doHadronScatter) ? "on" : "off") << endl
       << " afterDecay     = " << ((afterDecay)      ? "on" : "off") << endl
       << " allowDecayProd = " << ((allowDecayProd)  ? "on" : "off") << endl
       << " scatterRepeat  = " << ((scatterRepeat)   ? "on" : "off") << endl
       << " tile           = " << ((doTile) ? "on" : "off") << endl
       << "  yMin          = " << yMin << endl
       << "  yMax          = " << yMax << endl
       << "  ytMax         = " << ytMax << endl
       << "  ytSize        = " << ytSize << endl
       << "  ptMax         = " << ptMax << endl
       << "  ptSize        = " << ptSize << endl
       << endl
       << " hadronSelect   = " << hadronSelect << endl
       << "  N             = " << Npar << endl
       << "  k             = " << kPar << endl
       << "  p             = " << pPar << endl
       << endl
       << " scatterProb    = " << scatterProb << endl
       << "  j             = " << jPar << endl
       << "  rMax          = " << rMax << endl
       << endl
       << " pTsigma        = " << pTsigma2 << endl
       << " pT0MPI         = " << pT0MPI << endl
       << endl
       << " sigElMax       = " << sigElMax << endl << endl;

  return;
}


//--------------------------------------------------------------------------

// Perform hadron scattering

void HadronScatter::scatter(Event& event) {
  // Reset tiles
  for (int yt = 0; yt < ytMax; yt++)
    for (int pt = 0; pt < ptMax; pt++)
      tile[yt][pt].clear();

  // Generate list of hadrons which can take part in the scattering
  for (int i = 0; i < event.size(); i++)
    if (event[i].isFinal() && event[i].isHadron() && canScatter(event, i))
      tile[yTile(event, i)][pTile(event, i)].insert(HSIndex(i, i));

  // Generate all pairwise interaction probabilities
  vector < HadronScatterPair > scatterList;
  // For each tile and for each hadron in the tile do pairing
  for (int pt1 = 0; pt1 < ptMax; pt1++)
    for (int yt1 = 0; yt1 < ytMax; yt1++)
      for (set < HSIndex >::iterator si1 = tile[yt1][pt1].begin();
           si1 != tile[yt1][pt1].end(); si1++)
        tileIntProb(scatterList, event, *si1, yt1, pt1, false);
  // Sort by ordering measure (largest to smallest)
  sort(scatterList.rbegin(), scatterList.rend());

  // Reset list of things that have scattered
  if (scatterRepeat) scattered.clear();

  // Do scatterings
  while (scatterList.size() > 0) {
    // Check still valid and scatter
    HadronScatterPair &hsp = scatterList[0];
    if (!event[hsp.i1.second].isFinal() || !event[hsp.i2.second].isFinal()) {
      scatterList.erase(scatterList.begin());
      continue;
    }
    // Remove old entries in tiles and scatter
    tile[hsp.yt1][hsp.pt1].erase(hsp.i1);
    tile[hsp.yt2][hsp.pt2].erase(hsp.i2);
    hadronScatter(event, hsp);
    // Store new indices for later
    HSIndex iNew1 = hsp.i1, iNew2 = hsp.i2;

    // Check if hadrons can scatter again
    bool resort = false;
    if (scatterRepeat) {
      // Check for new scatters of iNew1 and iNew2
      HSIndex iNew = iNew1;
      for (int i = 0; i < 2; i++) {
        if (canScatter(event, iNew.second)) {
          // If both can scatter again, make sure they can't scatter
          // with each other
          if (i == 1) scattered.insert(HSIndex(min(iNew1.first, iNew2.first),
                                             max(iNew1.first, iNew2.first)));
          int yt = yTile(event, iNew.second);
          int pt = pTile(event, iNew.second);
          resort = tileIntProb(scatterList, event, iNew, yt, pt, true);
          tile[yt][pt].insert(iNew);
        }
        iNew = iNew2;
      } // for (i)

    } // if (scatterRepeat)

    // Remove the old entry and onto next
    scatterList.erase(scatterList.begin());
    // Resort list if anything added
    if (resort) sort(scatterList.rbegin(), scatterList.rend());

  } // while (scatterList.size() > 0)

  // Done.
  return;
}


//--------------------------------------------------------------------------

// Criteria for if a hadron is allowed to scatter or not

bool HadronScatter::canScatter(Event& event, int i) {
  // Probability to accept
  double p = 0.;

  // Pions, K+, K-, p+, pbar- only
  if (scatterProb == 1 || scatterProb == 2)
    if (event[i].idAbs() != 111 && event[i].idAbs() != 211 &&
        event[i].idAbs() != 321 && event[i].idAbs() != 2212)
      return false;

  // Probability
  switch (hadronSelect) {
  case 0:
    double t1 = exp( - event[i].pT2() / 2. / pTsigma2);
    double t2 = pow(pT0MPI, pPar) /
                pow(pT0MPI * pT0MPI + event[i].pT2(), pPar / 2.);
    p = Npar * t1 / ( (1 - kPar) * t1 + kPar * t2 );
    break;
  }

  // Return true/false
  if (p > rndmPtr->flat()) return true;
  else                     return false;
}


//--------------------------------------------------------------------------

// Probability for scattering
bool HadronScatter::doesScatter(Event& event, const HSIndex &i1,
  const HSIndex &i2) {
  Particle &p1 = event[i1.second];
  Particle &p2 = event[i2.second];

  // Can't come from the same decay
  if (!allowDecayProd && event[i1.first].mother1() ==
      event[i2.first].mother1() && event[event[i1.first].mother1()].isHadron())
    return false;

  // Check that the two hadrons have not already scattered
  if (scatterRepeat &&
      scattered.find(HSIndex(min(i1.first, i2.first), max(i1.first, i2.first)))
      != scattered.end()) return false;

  // K-K, p-p and K-p not allowed
  int id1 = min(p1.idAbs(), p2.idAbs());
  int id2 = max(p1.idAbs(), p2.idAbs());
  if (scatterProb == 1 || scatterProb == 2) {
    if ((id1 == 321 || id1 == 2212) && id1 == id2) return false;
    if (id1 == 321 && id2 == 2212)                 return false;
  }

  // Distance in y - phi space
  double dy  = p1.y() - p2.y();
  double dp  = abs(p1.phi() - p2.phi());
  if (dp > M_PI) dp = 2 * M_PI - dp;
  double dr2 = dy * dy + dp * dp;
  double p   = max(0., 1. - dr2 / rMax2);

  // Additional factor depending on scatterProb
  if (scatterProb == 0 || scatterProb == 1) {
    p *= jPar;

  // Additional pair dependent probability
  } else if (scatterProb == 2) {
    // Wcm and which partial wave amplitude to use
    double Wcm = (p1.p() + p2.p()).mCalc();
    int    pw = 0;
    if      ((id1 == 111 || id1 == 211) && (id2 == 111 || id2 == 211)) pw = 0;
    else if ((id1 == 111 || id1 == 211) && id2 == 321)                 pw = 1;
    else if ((id1 == 111 || id1 == 211) && id2 == 2212)                pw = 2;
    else {
      infoPtr->errorMsg("Error in HadronScatter::doesScatter:"
        "unknown subprocess");
    }
    if (!sigmaPW[pw].setSubprocess(p1.id(), p2.id())) {
      infoPtr->errorMsg("Error in HadronScatter::doesScatter:"
        "setSubprocess failed");
    } else {
      p *= 1 - exp(-jPar * sigmaPW[pw].sigmaEl(Wcm));
    }
  }

  // Return true/false
  return (p > rndmPtr->flat());
}


//--------------------------------------------------------------------------

// Calculate ordering measure

double HadronScatter::measure(Event& event, int idx1, int idx2) {
  Particle &p1 = event[idx1];
  Particle &p2 = event[idx2];
  return abs(p1.pT() / p1.mT() - p2.pT() / p2.mT());
}


//--------------------------------------------------------------------------

// Scatter a pair

bool HadronScatter::hadronScatter(Event& event, HadronScatterPair &hsp) {
  bool doSwap = (0.5 < rndmPtr->flat()) ? true : false;
  if (doSwap) swap(hsp.i1, hsp.i2);

  Particle &p1 = event[hsp.i1.second];
  Particle &p2 = event[hsp.i2.second];

  // Pick theta and phi (always flat)
  double ct = 0., phi = 2 * M_PI * rndmPtr->flat();

  // Flat for all flavours
  if (scatterProb == 0 || scatterProb == 1) {
    ct  = 2. * rndmPtr->flat() - 1.;

  // pi-pi, pi-K and pi-p only using partial wave amplitudes
  } else if (scatterProb == 2) {
    // Wcm and which partial wave amplitude to use
    int    id1 = min(p1.idAbs(), p2.idAbs());
    int    id2 = max(p1.idAbs(), p2.idAbs());
    double Wcm = (p1.p() + p2.p()).mCalc();
    int    pw  = 0;
    if      ((id1 == 111 || id1 == 211) && (id2 == 111 || id2 == 211)) pw = 0;
    else if ((id1 == 111 || id1 == 211) && id2 == 321)                 pw = 1;
    else if ((id1 == 111 || id1 == 211) && id2 == 2212)                pw = 2;
    else {
      infoPtr->errorMsg("Error in HadronScatter::hadronScatter:"
        "unknown subprocess");
    }
    sigmaPW[pw].setSubprocess(p1.id(), p2.id());
    ct = sigmaPW[pw].pickCosTheta(Wcm);
  }

  // Rotation
  RotBstMatrix sMat;
  sMat.toCMframe(p1.p(), p2.p());
  sMat.rot(acos(ct), phi);
  sMat.fromCMframe(p1.p(), p2.p());
  Vec4 v1 = p1.p(), v2 = p2.p();
  v1.rotbst(sMat);
  v2.rotbst(sMat);

  // Update event record
  int iNew1 = event.copy(hsp.i1.second, 98);
  event[iNew1].p(v1);
  event[iNew1].e(event[iNew1].eCalc());
  event[hsp.i1.second].statusNeg();
  int iNew2 = event.copy(hsp.i2.second, 98);
  event[iNew2].p(v2);
  event[iNew2].e(event[iNew2].eCalc());
  event[hsp.i2.second].statusNeg();

  // New indices in event record
  hsp.i1.second = iNew1;
  hsp.i2.second = iNew2;

  if (doSwap) swap(hsp.i1, hsp.i2);
  return true;
}


//--------------------------------------------------------------------------

// Calculate pair interaction probabilities in nearest-neighbour tiles
// (yt1, pt1) represents centre cell (8):
//    7 | 0 | 1
//   -----------
//    6 | 8 | 2
//   -----------
//    5 | 4 | 3

bool HadronScatter::tileIntProb(vector < HadronScatterPair > &hsp,
    Event &event, const HSIndex &i1, int yt1, int pt1, bool doAll) {
  // Track if a new pair is added
  bool pairAdded = false;

  // Special handling for pairing in own tile
  if (!doAll) {
    set < HSIndex >::iterator si2 = tile[yt1][pt1].find(i1);
    while (++si2 != tile[yt1][pt1].end())
      // Calculate if scatters
      if (doesScatter(event, i1, *si2)) {
        double m = measure(event, i1.second, si2->second);
        hsp.push_back(HadronScatterPair(i1, yt1, pt1, *si2, yt1, pt1, m));
      }
  }

  // And the rest
  int tileMax = (doAll) ? 9 : 4;
  for (int tileNow = 0; tileNow < tileMax; tileNow++) {
    // New (yt, pt) coordinate
    int yt2 = yt1, pt2 = pt1;
    switch (tileNow) {
    case 0:        ++pt2; break;
    case 1: ++yt2; ++pt2; break;
    case 2: ++yt2;        break;
    case 3: ++yt2; --pt2; break;
    case 4:        --pt2; break;
    case 5: --yt2; --pt2; break;
    case 6: --yt2;        break;
    case 7: --yt2; ++pt2; break;
    }

    // Limit in rapidity tiles
    if (yt2 < 0 || yt2 >= ytMax) continue;
    // Wraparound for phi tiles
    if (pt2 < 0) {
      pt2 = ptMax - 1;
      if (pt2 == pt1 || pt2 == pt1 + 1) continue;

    } else if (pt2 >= ptMax) {
      pt2 = 0;
      if (pt2 == pt1 || pt2 == pt1 - 1) continue;
    }

    // Interaction probability
    for (set < HSIndex >::iterator si2 = tile[yt2][pt2].begin();
         si2 != tile[yt2][pt2].end(); si2++) {
      // Calculate if scatters
      if (doesScatter(event, i1, *si2)) {
        double m = measure(event, i1.second, si2->second);
        hsp.push_back(HadronScatterPair(i1, yt1, pt1, *si2, yt2, pt2, m));
        pairAdded = true;
      }
    }
  }

  return pairAdded;
}

//==========================================================================

} // namespace Pythia8
