// DeuteronProduction.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// DeuteronProduction class.

#include "Pythia8/DeuteronProduction.h"

namespace Pythia8 {

//==========================================================================

// The DeuteronProduction class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of times to try a decay sampling.
const int DeuteronProduction::NTRYDECAY = 10;

  // These numbers are hardwired empirical parameters,
// intended to speed up the M-generator.
const double DeuteronProduction::WTCORRECTION[11] = { 1., 1., 1.,
  2., 5., 15., 60., 250., 1250., 7000., 50000. };

//--------------------------------------------------------------------------

// Find settings. Precalculate table used to find momentum shifts.

bool DeuteronProduction::init(Info* infoPtrIn, Settings& settings,
 ParticleData* pdbPtrIn, Rndm* rndmPtrIn) {

  // Save the pointers.
  infoPtr = infoPtrIn;
  pdbPtr  = pdbPtrIn;
  rndmPtr = rndmPtrIn;

  // Parse the settings.
  valid = true; ids.clear(); parms.clear(); masses.clear();
  models = settings.mvec("DeuteronProduction:models");
  vector<string> wvec;
  wvec = settings.wvec("DeuteronProduction:channels");
  for (int iWvec = 0; iWvec < int(wvec.size()); ++iWvec)
    ids.push_back(parseIds(wvec[iWvec]));
  wvec = settings.wvec("DeuteronProduction:parms");
  for (int iWvec = 0; iWvec < int(wvec.size()); ++iWvec)
    parms.push_back(parseParms(wvec[iWvec]));

  // Technical settings.
  bool verbose(settings.flag("Init:showProcesses"));
  mSafety = settings.parm("ParticleDecays:mSafety");
  kMin    = settings.parm("DeuteronProduction:kMin");
  kMax    = settings.parm("DeuteronProduction:kMax");
  kTol    = settings.parm("DeuteronProduction:kTol");
  kSteps  = settings.mode("DeuteronProduction:kSteps");

  // Check the configuration vectors.
  string pre("Error in DeuteronProduction::init: ");
  if (parms.size() != ids.size() || parms.size() != models.size()) {
    ostringstream vals;
    vals << " " << ids.size() << ", " << models.size() << ", " << parms.size();
    infoPtr->errorMsg(pre + "channels, models, and parms are size",
      vals.str());
    valid = false;
  }

  // Check each channel configuration.
  for (int chn = 0; chn < int(parms.size()); ++chn) {
    if (ids[chn].size() < 3) {
      infoPtr->errorMsg(pre + "ids must have 3 or more IDs");
      valid = false;
    } else if (ids[chn][2] != 0) {
      infoPtr->errorMsg(pre + "ids initial state must be size 2");
      valid = false;
    } if (models[chn] == 0 && parms[chn].size() != 2) {
      infoPtr->errorMsg(pre + "model 0 channels must have", "2 coefficients");
      valid = false;
    } if (models[chn] == 1 && parms[chn].size() != 15) {
      infoPtr->errorMsg(pre + "model 1 channels must have", "15 coefficients");
      valid = false;
    } if (models[chn] == 2 && parms[chn].size() != 5) {
      infoPtr->errorMsg(pre + "model 2 channels must have", "2 coefficients");
      valid = false;
    } if (models[chn] == 3 && parms[chn].size()%5 != 0) {
      infoPtr->errorMsg(pre + "model 3 channels must have",
        "a multiple of 5 coefficients");
      valid = false;
    }
  }
  if (!valid) return valid;
  mPion = pdbPtr->m0(211);

  // Find channel maxima and set the normalization.
  if (verbose)
    cout << "\n *----------  PYTHIA Deuteron Production "
         << "Initialization  -----------*\n"
         << " |" << setw(68) << "|           |           |\n"
         << " | Subprocess" << setw(57) << "| k (GeV)   | Max (mb)  |\n"
         << " |" << setw(68) << "|           |           |\n"
         << " |--------------------------------------------------------------"
         << "----|\n"
         << " |" << setw(68) << "|           |           |\n";
  double max(0), k, s;
  for (int chn = 0; chn < int(ids.size()); ++chn) {

    // Always require the proton first and add the deuteron.
    if (ids[chn][1] == 2212) swap(ids[chn][1], ids[chn][0]);
    ids[chn].push_back(1000010020);

    // Set the nominal masses.
    vector<double> mass(ids[chn].size(), 0);
    for (int id = 0; id < int(ids[chn].size()); ++id)
      mass[id] = pdbPtr->m0(ids[chn][id]);
    masses.push_back(mass);

    // Calculate the maximum cross-section.
    maximum(k, s, chn);
    if (verbose) {
      string proc(" |");
      for (int id = 0; id < 2; ++id) proc += " " + pdbPtr->name(ids[chn][id]);
      proc += " ->";
      for (int id = 3; id < int(ids[chn].size()); ++id)
        proc += " " + pdbPtr->name(ids[chn][id]);
      cout << left << setw(43) << proc;
      cout << " | " << scientific << setprecision(3) << k
           << " | " << scientific << setprecision(3) << s << " |\n";
    }
    if (s > max) max = s;
  }

  // Set normalization.
  norm = settings.parm("DeuteronProduction:norm");
  if (norm < 1) norm = max;
  else norm *= max;
  if (verbose)
    cout << " |" << right << setw(68) << "                        |\n"
         << " | Using a maximum of " << norm << " mb" << setw(36) << "|\n"
         << " |" << right << setw(68) << "                        |\n"
         << " *----------  End PYTHIA Deuteron Production "
         << "Initializaiton  -------*\n";

  // Done.
  return valid;

}

//--------------------------------------------------------------------------

// Form deuterons in an event.

bool DeuteronProduction::combine(Event& event) {

  // Create nucleon and anti-nucleon vectors.
  if (!valid) return false;
  vector<int> nucs, anucs;
  for (int idx = 0; idx < event.size(); ++idx) {
    Particle &prt = event[idx];
    if (prt.statusAbs() > 80 && (prt.idAbs() == 2212 || prt.idAbs() == 2112)
        && prt.iBotCopy() == idx) {
      if (prt.id() > 0) nucs.push_back(idx);
      else anucs.push_back(idx);
      prt.undoDecay();
    }
  }

  // Bind the combinations, make used nucleon energies positive, and return.
  bind(event, nucs);
  bind(event, anucs);
  return true;
}

//--------------------------------------------------------------------------

// Bind the nucleon-pair combinations.

void DeuteronProduction::bind(Event& event, vector<int>& prts) {

  // Build the combinations and set the cross-section for each channel.
  vector<pair<int, int> > cmbs;
  combos(event, prts, cmbs);
  vector<double> sigmas(ids.size(), 0);

  // Loop over the nucleon pair combinations.
  for (int cmb = 0; cmb < int(cmbs.size()); ++ cmb) {

    // Skip if the pair has already been bound.
    Particle &prt0 = event[cmbs[cmb].first];
    Particle &prt1 = event[cmbs[cmb].second];
    if (prt0.status() < 0 || prt1.status() < 0) continue;

    // Calculate the momentum difference.
    Vec4 p0(prt0.p()), p1(prt1.p()), p(p0 + p1);
    p0.bstback(p);
    p1.bstback(p);
    double k((p0 - p1).pAbs());

    // Try binding each channel.
    double sum(0);
    for (int chn = 0; chn < int(ids.size()); ++chn) {
      if (prt0.idAbs() == ids[chn][0] && prt1.idAbs() == ids[chn][1])
        sigmas[chn] = sigma(k, chn);
      else {sigmas[chn] = 0; continue;}
      if (sigmas[chn] > norm)
        infoPtr->errorMsg("Warning in DeuteronProduction::bind:",
          "maximum weight exceeded");
      if (rndmPtr->flat() >= sigmas[chn]/norm) sigmas[chn] = 0;
      sum += sigmas[chn];
    }

    // Pick a bound channel.
    if (sum == 0) continue;
    double rndm(sum*rndmPtr->flat()); int chn(-1);
    do rndm -= sigmas[++chn];
    while (rndm > 0. && chn < int(sigmas.size()));
    if (chn < 0) continue;

    // Generate the decay and add to the event record.
    decay(event, prt0.index(), prt1.index(), chn);
  }

}

//--------------------------------------------------------------------------

// Build the nucleon-pair combinations and shuffle.

void DeuteronProduction::combos(Event& event, vector<int>& prts,
  vector<pair<int, int> >& cmbs) {

  // Create the combos.
  for (int idx0 = 0; idx0 < int(prts.size()); ++idx0) {
    int prt0(prts[idx0]), id(event[prt0].idAbs() == 2112);
    for (int idx1 = idx0 + 1; idx1 < int(prts.size()); ++idx1) {
      int prt1(prts[idx1]);
      cmbs.push_back(make_pair(id ? prt1 : prt0, id ? prt0 : prt1));
    }
  }

  // Shuffle.
  for (int idx = int(cmbs.size()) - 1; idx > 0; --idx)
    swap(cmbs[idx], cmbs[rndmPtr->flat()*(idx + 1)]);
}

//--------------------------------------------------------------------------

// Single pion final state fit, equations 10/13/14 of arXiv:1504.07242.

double DeuteronProduction::fit(double k, vector<double>& c, int i) {

  return c[i] * pow(k, c[i + 1])
    / (pow((c[i + 2] - exp(c[i + 3]*k)),2) + c[i + 4]);

}

//--------------------------------------------------------------------------

// Return the cross-section for a given channel.

double DeuteronProduction::sigma(double k, int chn) {

  double sum(0);
  int model(models[chn]);
  vector<double> &c = parms[chn];
  vector<double> &m = masses[chn];

  // Check allowed phase-space.
  double ecm(sqrt(m[0]*m[0] + k*k/4) + sqrt(m[1]*m[1] + k*k/4)), mtot(0);
  for (int dtr = 3; dtr < int(m.size()); ++dtr) mtot += m[dtr];
  if (ecm < mtot) {
    sum = 0;

  // Step function, e.g. coalescence model.
  } else if (model == 0) {
    sum = k < c[0] ? c[1] : 0;

  // p n -> gamma d, equation 7 where first parameter is function k-split.
  } else if (model == 1) {
    if (k < c[0]) {for (int i = 1; i < 13; ++i) sum += c[i]*pow(k, i - 2);}
    else sum = exp(-c[13]*k - c[14]*k*k);

  // p/n p/n -> pi d, equation 10.
  } else if (model == 2)  {
    double s(ecm*ecm), q(sqrtpos(pow(s + m[3]*m[3] - m.back()*m.back(), 2)
                                 /(4*s) - m[3]*m[3]));
    sum = fit(q/mPion, c, 0);

  // p/n p/n -> pi pi d, equations 13 and 14.
  } else if (model == 3) {
    for (int i = 0; i < int(c.size()); i += 5) sum += fit(k, c, i);
  }
  return sum*1e-3;

}

//--------------------------------------------------------------------------

// N-body decay using the M-generator algorithm described in "Monte
// Carlo Phase Space" by F. James in CERN 68-15, May 1968. Modified
// from ParticleDecays::mGenerator.

bool DeuteronProduction::decay(Event& event, int idx0, int idx1, int chn) {

  // Set the decay product masses and fill the event.
  if (idx0 < idx1) swap(idx0, idx1);
  int mult(ids[chn].size() - 3);
  vector<double> mProd(mult + 1), mInv(mult + 1, 0);
  Vec4 pMom(event[idx0].p() + event[idx1].p());
  mProd[0] = pMom.mCalc();

  // Select masses. Fail if too close or inconsistent.
  double mDiff(0);
  for (int tries = 0; tries < NTRYDECAY && mDiff < mSafety; ++tries) {
    mDiff = mProd[0];
    for (int i = 1; i <= mult; ++i) {
      mProd[i] = pdbPtr->mSel(ids[chn][i + 2]);
      mDiff -= mProd[i];
    }
  }
  if (mDiff < mSafety) {
    infoPtr->errorMsg("Warning in DeuteronProduction::decay:",
                      "no valid decay found");
    return false;
  }

  // Set up the event and invariant masses.
  vector<int> iProd(mult + 1);
  for (int i = 1; i <= mult; ++i) {
    int id(ids[chn][i + 2]);
    if (event[idx0].id() < 0 && pdbPtr->hasAnti(id)) id *= -1;
    iProd[i] = event.append(id, 121, idx0, idx1, 0, 0, 0, 0,
                            Vec4(0., 0., 0., 0.), mProd[i], 0);
  }
  for (int i = 0; i <= mult; ++i) mInv[i] = mProd[i];
  vector<double> rndmOrd(mult, 0);
  vector<Vec4> pInv(mult + 1, 0);

  // Calculate the maximum weight in the decay.
  double wtPS, wtME, wtMEmax;
  double wtPSmax = 1. / WTCORRECTION[mult];
  double mMax    = mDiff + mProd[mult];
  double mMin    = 0.;
  for (int i = mult - 1; i > 0; --i) {
    mMax        += mProd[i];
    mMin        += mProd[i+1];
    double mNow  = mProd[i];
    wtPSmax     *= 0.5 * sqrtpos( (mMax - mMin - mNow) * (mMax + mMin + mNow)
                 * (mMax + mMin - mNow) * (mMax - mMin + mNow) ) / mMax;
  }

  // Begin loop over matrix-element corrections.
  do {
    wtME    = 1.;
    wtMEmax = 1.;

    // Begin loop to find the set of intermediate invariant masses.
    do {
      wtPS  = 1.;

      // Find and order random numbers in descending order.
      rndmOrd[0] = 1.;
      for (int i = 1; i < mult - 1; ++i) {
        double rndm = rndmPtr->flat();
        rndmOrd[i] = rndm;
        for (int j = i - 1; j > 0; --j) {
          if (rndm > rndmOrd[j]) swap( rndmOrd[j], rndmOrd[j+1] );
          else break;
        }
      }
      rndmOrd[mult -1] = 0.;

      // Translate into intermediate masses and find weight.
      for (int i = mult - 1; i > 0; --i) {
        mInv[i] = mInv[i+1] + mProd[i] + (rndmOrd[i-1] - rndmOrd[i]) * mDiff;
        wtPS   *= 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - mProd[i])
          * (mInv[i] + mInv[i+1] + mProd[i]) * (mInv[i] + mInv[i+1] - mProd[i])
          * (mInv[i] - mInv[i+1] + mProd[i]) ) / mInv[i];
      }

    // If rejected, try again with new invariant masses.
    } while ( wtPS < rndmPtr->flat() * wtPSmax );

    // Perform two-particle decays in the respective rest frame.
    for (int i = 1; i < mult; ++i) {
      double pAbs = 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - mProd[i])
        * (mInv[i] + mInv[i+1] + mProd[i]) * (mInv[i] + mInv[i+1] - mProd[i])
        * (mInv[i] - mInv[i+1] + mProd[i]) ) / mInv[i];

      // Isotropic angles give three-momentum.
      double cosTheta = 2. * rndmPtr->flat() - 1.;
      double sinTheta = sqrt(1. - cosTheta*cosTheta);
      double phi      = 2. * M_PI * rndmPtr->flat();
      double pX       = pAbs * sinTheta * cos(phi);
      double pY       = pAbs * sinTheta * sin(phi);
      double pZ       = pAbs * cosTheta;

      // Calculate energies, fill four-momenta.
      double eHad     = sqrt( mProd[i]*mProd[i] + pAbs*pAbs);
      double eInv     = sqrt( mInv[i+1]*mInv[i+1] + pAbs*pAbs);
      event[iProd[i]].p( pX, pY, pZ, eHad);
      pInv[i+1].p( -pX, -pY, -pZ, eInv);
    }

    // Boost decay products to the mother rest frame.
    event[iProd[mult]].p( pInv[mult] );
    for (int iFrame = mult - 1; iFrame > 1; --iFrame)
      for (int i = iFrame; i <= mult; ++i)
        event[iProd[i]].bst( pInv[iFrame], mInv[iFrame]);

  // If rejected, try again with new invariant masses.
  } while ( wtME < rndmPtr->flat() * wtMEmax );

  // Boost decay products to the current frame and set nucleon properties.
  for (int i = 1; i <= mult; ++i) event[iProd[i]].bst( pMom, mInv[1] );
  event[idx0].statusNeg();
  event[idx1].statusNeg();
  event[idx0].daughter1(iProd[1]); event[idx0].daughter2(iProd.back());
  event[idx1].daughter1(iProd[1]); event[idx1].daughter2(iProd.back());
  return true;

}

//--------------------------------------------------------------------------

// Find the maximum for a cross-section.

void DeuteronProduction::maximum(double& k, double& s, int chn) {

  // Perform a grid search.
  double y, xa(kMin), xb(kMax), xstep((xb - xa)/(kSteps + 1)), xm(xa), ym(0);
  for (double x = xa; x <= xb; x += xstep) {
    y = sigma(x, chn);
    if (y > ym) {xm = x; ym = y;}
  }

  // Run the simplex algorithm.
  vector<double> xs(5, xm); int im(2);
  xs[0] = xm == xa ? xa : xm - xstep;
  xs[4] = xm == xb ? xb : xm + xstep;
  for (int itr = 0; itr < 1000 && abs((xs[0] - xs[4])/xs[2]) > kTol; ++itr) {
    xs[2] = (xs[0] + xs[4])/2;
    xs[1] = (xs[0] + xs[2])/2;
    xs[3] = (xs[2] + xs[4])/2;
    im = 0;
    for (int i = 0; i < int(xs.size()); ++i) {
      y = sigma(xs[i], chn);
      if (y > ym) {im = i; ym = y;}
    }
    if (im < 2) xs[4] = xs[2];
    else if (im > 2) xs[0] = xs[2];
    else {xs[0] = xs[1]; xs[4] = xs[3];}
  }
  k = xs[im]; s = ym;

}

//--------------------------------------------------------------------------

// Parse the IDs.

vector<int> DeuteronProduction::parseIds(string line) {

  vector<int> vals;
  if (line == "") return vals;
  int val;
  size_t pos(0);
  while (pos != string::npos) {
    pos = line.find(" ");
    if (pos == 0) {line = line.substr(pos + 1); continue;}
    istringstream stream(line.substr(0, pos));
    line = line.substr(pos + 1);
    stream >> val;
    vals.push_back(val);
  }
  return vals;

}

//--------------------------------------------------------------------------

// Parse the parameters.

vector<double> DeuteronProduction::parseParms(string line) {

  vector<double> vals;
  if (line == "") return vals;
  double val;
  size_t pos(0);
  while (pos != string::npos) {
    pos = line.find(" ");
    if (pos == 0) {line = line.substr(pos + 1); continue;}
    istringstream stream(line.substr(0, pos));
    line = line.substr(pos + 1);
    stream >> val;
    vals.push_back(val);
  }
  return vals;

}

//==========================================================================

} // end namespace Pythia8
