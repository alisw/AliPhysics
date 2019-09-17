// SigmaOnia.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// charmonia/bottomonia simulation classes.

#include "Pythia8/SigmaOnia.h"
#include <limits>
namespace Pythia8 {

//==========================================================================

// SigmaOniaSetup class.
// A helper class used to setup the SigmaOnia processes.

//--------------------------------------------------------------------------

// The constructor.

SigmaOniaSetup::SigmaOniaSetup(Info* infoPtrIn, Settings* settingsPtrIn,
  ParticleData* particleDataPtrIn, int flavourIn)
  : valid3S1(true), valid3PJ(true), valid3DJ(true), validDbl3S1(true),
  flavour(flavourIn) {

  // Set the pointers and category/key strings and mass splitting.
  infoPtr = infoPtrIn;
  settingsPtr = settingsPtrIn;
  particleDataPtr = particleDataPtrIn;
  cat   = (flavour == 4) ? "Charmonium" : "Bottomonium";
  key   = (flavour == 4) ? "ccbar" : "bbbar";
  mSplit = settingsPtr->parm("Onia:massSplit");
  if (!settingsPtr->flag("Onia:forceMassSplit")) mSplit = -mSplit;

  // Set the general switch settings.
  onia        = settingsPtr->flag("Onia:all");
  onia3S1     = settingsPtr->flag("Onia:all(3S1)");
  onia3PJ     = settingsPtr->flag("Onia:all(3PJ)");
  onia3DJ     = settingsPtr->flag("Onia:all(3DJ)");
  oniaFlavour = settingsPtr->flag(cat + ":all");

  // Set the names of the matrix element settings.
  meNames3S1.push_back(cat + ":O(3S1)[3S1(1)]");
  meNames3S1.push_back(cat + ":O(3S1)[3S1(8)]");
  meNames3S1.push_back(cat + ":O(3S1)[1S0(8)]");
  meNames3S1.push_back(cat + ":O(3S1)[3P0(8)]");
  meNames3PJ.push_back(cat + ":O(3PJ)[3P0(1)]");
  meNames3PJ.push_back(cat + ":O(3PJ)[3S1(8)]");
  meNames3DJ.push_back(cat + ":O(3DJ)[3D1(1)]");
  meNames3DJ.push_back(cat + ":O(3DJ)[3P0(8)]");
  meNamesDbl3S1.push_back(cat + ":O(3S1)[3S1(1)]1");
  meNamesDbl3S1.push_back(cat + ":O(3S1)[3S1(1)]2");

  // Set the names of the production settings.
  ggNames3S1.push_back(cat + ":gg2" + key + "(3S1)[3S1(1)]g");
  ggNames3S1.push_back(cat + ":gg2" + key + "(3S1)[3S1(1)]gm");
  ggNames3S1.push_back(cat + ":gg2" + key + "(3S1)[3S1(8)]g");
  ggNames3S1.push_back(cat + ":gg2" + key + "(3S1)[1S0(8)]g");
  ggNames3S1.push_back(cat + ":gg2" + key + "(3S1)[3PJ(8)]g");
  qgNames3S1.push_back(cat + ":qg2" + key + "(3S1)[3S1(8)]q");
  qgNames3S1.push_back(cat + ":qg2" + key + "(3S1)[1S0(8)]q");
  qgNames3S1.push_back(cat + ":qg2" + key + "(3S1)[3PJ(8)]q");
  qqNames3S1.push_back(cat + ":qqbar2" + key + "(3S1)[3S1(8)]g");
  qqNames3S1.push_back(cat + ":qqbar2" + key + "(3S1)[1S0(8)]g");
  qqNames3S1.push_back(cat + ":qqbar2" + key + "(3S1)[3PJ(8)]g");
  ggNames3PJ.push_back(cat + ":gg2" + key + "(3PJ)[3PJ(1)]g");
  ggNames3PJ.push_back(cat + ":gg2" + key + "(3PJ)[3S1(8)]g");
  qgNames3PJ.push_back(cat + ":qg2" + key + "(3PJ)[3PJ(1)]q");
  qgNames3PJ.push_back(cat + ":qg2" + key + "(3PJ)[3S1(8)]q");
  qqNames3PJ.push_back(cat + ":qqbar2" + key + "(3PJ)[3PJ(1)]g");
  qqNames3PJ.push_back(cat + ":qqbar2" + key + "(3PJ)[3S1(8)]g");
  ggNames3DJ.push_back(cat + ":gg2" + key + "(3DJ)[3DJ(1)]g");
  ggNames3DJ.push_back(cat + ":gg2" + key + "(3DJ)[3PJ(8)]g");
  qgNames3DJ.push_back(cat + ":qg2" + key + "(3DJ)[3PJ(8)]q");
  qqNames3DJ.push_back(cat + ":qqbar2" + key + "(3DJ)[3PJ(8)]g");
  dblNames3S1.push_back(cat + ":gg2double" + key + "(3S1)[3S1(1)]");
  dblNames3S1.push_back(cat + ":qqbar2double" + key + "(3S1)[3S1(1)]");

  // Initialise and check all settings.
  states3S1 = settingsPtr->mvec(cat + ":states(3S1)");
  initStates("(3S1)", states3S1, spins3S1, valid3S1);
  initSettings("(3S1)", states3S1.size(), meNames3S1, mes3S1, valid3S1);
  initSettings("(3S1)", states3S1.size(), ggNames3S1, ggs3S1, valid3S1);
  initSettings("(3S1)", states3S1.size(), qgNames3S1, qgs3S1, valid3S1);
  initSettings("(3S1)", states3S1.size(), qqNames3S1, qqs3S1, valid3S1);
  states3PJ = settingsPtr->mvec(cat + ":states(3PJ)");
  initStates("(3PJ)", states3PJ, spins3PJ, valid3PJ);
  initSettings("(3PJ)", states3PJ.size(), meNames3PJ, mes3PJ, valid3PJ);
  initSettings("(3PJ)", states3PJ.size(), ggNames3PJ, ggs3PJ, valid3PJ);
  initSettings("(3PJ)", states3PJ.size(), qgNames3PJ, qgs3PJ, valid3PJ);
  initSettings("(3PJ)", states3PJ.size(), qqNames3PJ, qqs3PJ, valid3PJ);
  states3DJ = settingsPtr->mvec(cat + ":states(3DJ)");
  initStates("(3DJ)", states3DJ, spins3DJ, valid3DJ);
  initSettings("(3DJ)", states3DJ.size(), meNames3DJ, mes3DJ, valid3DJ);
  initSettings("(3DJ)", states3DJ.size(), ggNames3DJ, ggs3DJ, valid3DJ);
  initSettings("(3DJ)", states3DJ.size(), qgNames3DJ, qgs3DJ, valid3DJ);
  initSettings("(3DJ)", states3DJ.size(), qqNames3DJ, qqs3DJ, valid3DJ);
  states1Dbl3S1 = settingsPtr->mvec(cat + ":states(3S1)1");
  states2Dbl3S1 = settingsPtr->mvec(cat + ":states(3S1)2");
  initStates("(3S1)1", states1Dbl3S1, spins1Dbl3S1, validDbl3S1, false);
  initStates("(3S1)2", states2Dbl3S1, spins2Dbl3S1, validDbl3S1, false);
  if (states1Dbl3S1.size() != states2Dbl3S1.size()) {
    infoPtr->errorMsg("Error in SigmaOniaSetup: mvecs Charmonium:states(3S1)"
      " 1 and 2 are not the same size");
    validDbl3S1 = false;
    return;
  }
  initSettings("(3S1)1", states1Dbl3S1.size(), meNamesDbl3S1, mesDbl3S1,
    validDbl3S1);
  initSettings("(3S1)1", states1Dbl3S1.size(), dblNames3S1, dbls3S1,
    validDbl3S1);

}

//--------------------------------------------------------------------------

// Initialise the SigmaProcesses for g g -> X g production.

void SigmaOniaSetup::setupSigma2gg(vector<SigmaProcess*> &procs, bool oniaIn) {

  // Initialise the 3S1 processes.
  if (valid3S1) {
    for (unsigned int i = 0; i < states3S1.size(); ++i) {
      bool flag = oniaIn || onia || onia3S1 || oniaFlavour;
      // Colour-singlet.
      if (flag || ggs3S1[0][i])
        procs.push_back(new Sigma2gg2QQbar3S11g
          (states3S1[i], mes3S1[0][i], flavour*100 + 1));
      if (flag || ggs3S1[1][i])
        procs.push_back(new Sigma2gg2QQbar3S11gm
          (states3S1[i], mes3S1[0][i], flavour*110 + 1));
      // Colour-octet.
      if (flag || ggs3S1[2][i])
        procs.push_back(new Sigma2gg2QQbarX8g
          (states3S1[i], mes3S1[1][i], 0, mSplit, flavour*100+2));
      if (flag || ggs3S1[3][i])
        procs.push_back(new Sigma2gg2QQbarX8g
          (states3S1[i], mes3S1[2][i], 1, mSplit, flavour*100+5));
      if (flag || ggs3S1[4][i])
        procs.push_back(new Sigma2gg2QQbarX8g
          (states3S1[i], mes3S1[3][i], 2, mSplit, flavour*100+8));
    }
  }

  // Initialise the 3PJ processes.
  if (valid3PJ) {
    for (unsigned int i = 0; i < states3PJ.size(); ++i) {
      bool flag = oniaIn || onia || onia3PJ || oniaFlavour;
      // Colour-singlet.
      if (flag || ggs3PJ[0][i]) {
        procs.push_back(new Sigma2gg2QQbar3PJ1g
          (states3PJ[i], mes3PJ[0][i], spins3PJ[i], flavour*100 + 11));
      }
      // Colour-octet.
      if (flag || ggs3PJ[1][i])
        procs.push_back(new Sigma2gg2QQbarX8g
          (states3PJ[i], mes3PJ[1][i], 0, mSplit, flavour*100+14));
    }
  }

  // Initialise the 3DJ processes.
  if (valid3DJ) {
    for (unsigned int i = 0; i < states3DJ.size(); ++i) {
      bool flag = oniaIn || onia || onia3DJ || oniaFlavour;
      // Colour-singlet.
      if (flag || ggs3DJ[0][i]) {
        procs.push_back(new Sigma2gg2QQbar3DJ1g
          (states3DJ[i], mes3DJ[0][i], spins3DJ[i], flavour*100 + 17));
      }
      // Colour-octet.
      if (flag || ggs3DJ[1][i]) {
        procs.push_back(new Sigma2gg2QQbarX8g
          (states3DJ[i], mes3DJ[1][i], 2, mSplit, flavour*100+18));
      }
    }
  }

}

//--------------------------------------------------------------------------

// Initialise the SigmaProcesses for q g -> X q production.

void SigmaOniaSetup::setupSigma2qg(vector<SigmaProcess*> &procs, bool oniaIn) {

  // Initialise the 3S1 processes.
  if (valid3S1) {
    for (unsigned int i = 0; i < states3S1.size(); ++i) {
      bool flag = oniaIn || onia || onia3S1 || oniaFlavour;
      // Colour-octet.
      if (flag || qgs3S1[0][i])
        procs.push_back(new Sigma2qg2QQbarX8q
          (states3S1[i], mes3S1[1][i], 0, mSplit, flavour*100+3));
      if (flag || qgs3S1[1][i])
        procs.push_back(new Sigma2qg2QQbarX8q
          (states3S1[i], mes3S1[2][i], 1, mSplit, flavour*100+6));
      if (flag || qgs3S1[2][i])
        procs.push_back(new Sigma2qg2QQbarX8q
          (states3S1[i], mes3S1[3][i], 2, mSplit, flavour*100+9));
    }
  }

  // Initialise the 3PJ processes.
  if (valid3PJ) {
    for (unsigned int i = 0; i < states3PJ.size(); ++i) {
      bool flag = oniaIn || onia || onia3PJ || oniaFlavour;
      // Colour-singlet.
      if (flag || qgs3PJ[0][i])
        procs.push_back(new Sigma2qg2QQbar3PJ1q
          (states3PJ[i], mes3PJ[0][i], spins3PJ[i], flavour*100 + 12));
      // Colour-octet.
      if (flag || qgs3PJ[1][i])
        procs.push_back(new Sigma2qg2QQbarX8q
          (states3PJ[i], mes3PJ[1][i], 0, mSplit, flavour*100+15));
    }
  }

  // Initialise the 3DJ processes.
  if (valid3DJ) {
    for (unsigned int i = 0; i < states3DJ.size(); ++i) {
      bool flag = oniaIn || onia || onia3DJ || oniaFlavour;
      // Colour-octet.
      if (flag || qgs3DJ[0][i])
        procs.push_back(new Sigma2qg2QQbarX8q
          (states3DJ[i], mes3DJ[1][i], 2, mSplit, flavour*100+19));
    }
  }

}

//--------------------------------------------------------------------------

// Initialise the SigmaProcesses for q qbar -> X g production.

void SigmaOniaSetup::setupSigma2qq(vector<SigmaProcess*> &procs, bool oniaIn) {

  // Initialise the 3S1 processes.
  if (valid3S1) {
    for (unsigned int i = 0; i < states3S1.size(); ++i) {
      bool flag = oniaIn || onia || onia3S1 || oniaFlavour;
      // Colour-octet.
      if (flag || qqs3S1[0][i])
        procs.push_back(new Sigma2qqbar2QQbarX8g
          (states3S1[i], mes3S1[1][i], 0, mSplit, flavour*100+4));
      if (flag || qqs3S1[1][i])
        procs.push_back(new Sigma2qqbar2QQbarX8g
          (states3S1[i], mes3S1[2][i], 1, mSplit, flavour*100+7));
      if (flag || qqs3S1[2][i])
        procs.push_back(new Sigma2qqbar2QQbarX8g
          (states3S1[i], mes3S1[3][i], 2, mSplit, flavour*100+10));
    }
  }

  // Initialise the 3PJ processes.
  if (valid3PJ) {
    for (unsigned int i = 0; i < states3PJ.size(); ++i) {
      bool flag = oniaIn || onia || onia3PJ || oniaFlavour;
      // Colour-singlet.
      if (flag || qqs3PJ[0][i])
        procs.push_back(new Sigma2qqbar2QQbar3PJ1g
          (states3PJ[i], mes3PJ[0][i], spins3PJ[i], flavour*100 + 13));
      // Colour-octet.
      if (flag || qqs3PJ[1][i])
        procs.push_back(new Sigma2qqbar2QQbarX8g
          (states3PJ[i], mes3PJ[1][i], 0, mSplit, flavour*100+16));
    }
  }

  // Initialise the 3DJ processes.
  if (valid3DJ) {
    for (unsigned int i = 0; i < states3DJ.size(); ++i) {
      bool flag = oniaIn || onia || onia3DJ || oniaFlavour;
      // Colour-octet.
      if (flag || qqs3DJ[0][i])
        procs.push_back(new Sigma2qqbar2QQbarX8g
          (states3DJ[i], mes3DJ[1][i], 2, mSplit, flavour*100+20));
    }
  }

}

//--------------------------------------------------------------------------

// Initialise the SigmaProcesses for double onium production.

void SigmaOniaSetup::setupSigma2dbl(vector<SigmaProcess*> &procs,
  bool oniaIn) {

  // Initialise the 3S1 processes.
  if (validDbl3S1) {
    for (unsigned int i = 0; i < states1Dbl3S1.size(); ++i) {
      bool flag = oniaIn || onia || onia3S1 || oniaFlavour;
      // Colour-singlet.
      if ((flag || dbls3S1[0][i])) procs.push_back(
        new Sigma2gg2QQbar3S11QQbar3S11( states1Dbl3S1[i], states2Dbl3S1[i],
           mesDbl3S1[0][i], mesDbl3S1[1][i], flavour*100 + 21) );
      if ((flag || dbls3S1[1][i])) procs.push_back(
        new Sigma2qqbar2QQbar3S11QQbar3S11( states1Dbl3S1[i], states2Dbl3S1[i],
           mesDbl3S1[0][i], mesDbl3S1[1][i], flavour*100 + 22) );
    }
  }

}

//--------------------------------------------------------------------------

// Initialise and check the flavour, j-number, and validity of states.

void SigmaOniaSetup::initStates(string wave, const vector<int> &states,
  vector<int> &jnums, bool &valid, bool duplicates) {

  set<int> unique;
  unsigned int nstates(0);
  for (unsigned int i = 0; i < states.size(); ++i) {

    // Check state is unique and remove if not.
    stringstream state;
    state << states[i];
    unique.insert(states[i]);
    if (duplicates && nstates + 1 != unique.size()) {
      infoPtr->errorMsg("Error in SigmaOniaSetup::initStates: particle "
                        + state.str() + " in mvec " + cat + ":states"
                        + wave, "has duplicates");
      valid = false;
    } else ++nstates;

    // Determine quark composition and quantum numbers.
    int mod1(10), mod2(1);
    vector<int> digits;
    while (digits.size() < 7) {
      digits.push_back((states[i]%mod1 - states[i]%mod2) / mod2);
      mod1 *= 10;
      mod2 *= 10;
    }
    int s, l, j((digits[0] - 1)/2);
    if (j != 0) {
      if      (digits[4] == 0) {l = j - 1; s = 1;}
      else if (digits[4] == 1) {l = j;     s = 0;}
      else if (digits[4] == 2) {l = j;     s = 1;}
      else                     {l = j + 1; s = 1;}
    } else {
      if      (digits[4] == 0) {l = 0;  s = 0;}
      else                     {l = 1;  s = 1;}
    }

    // Check state validity.
    if (states[i] != 0) {
      if (!particleDataPtr->isParticle(states[i])) {
        infoPtr->errorMsg("Error in SigmaOniaSetup::initStates: particle "
                          + state.str() + " in mvec " + cat + ":states"
                          + wave, "is unknown");
        valid = false;
      }
      if (digits[3] != 0) {
        infoPtr->errorMsg("Error in SigmaOniaSetup::initStates: particle "
                          + state.str() + " in mvec " + cat + ":states"
                          + wave, " is not a meson");
        valid = false;
      }
      if (digits[2] != digits[1] || digits[1] != flavour) {
        infoPtr->errorMsg("Error in SigmaOniaSetup::initStates: particle "
                          + state.str() + " in mvec " + cat + ":states"
                          + wave, "is not a " + key + " state");
        valid = false;
      }
      if ((wave == "3S1" && (s != 1 || l != 0 || j != 1)) ||
          (wave == "3PJ" && (s != 1 || l != 1 || j < 0 || j > 2)) ||
          (wave == "3DJ" && (s != 1 || l != 2 || j < 1 || j > 3))) {
        infoPtr->errorMsg("Error in SigmaOniaSetup::initStates: particle "
                          + state.str() + " in mvec " + cat + ":states"
                          + wave, "is not a " + wave + " state");
        valid = false;
      }
    } else valid = false;
    jnums.push_back(j);
  }

}

//--------------------------------------------------------------------------

// Initialise and check a group of PVec settings.

void SigmaOniaSetup::initSettings(string wave, unsigned int size,
  const vector<string> &names, vector< vector<double> > &pvecs,
  bool &valid) {

  for (unsigned int i = 0; i < names.size(); ++i) {
    pvecs.push_back(settingsPtr->pvec(names[i]));
    if (pvecs.back().size() != size) {
      infoPtr->errorMsg("Error in SigmaOniaSetup::initSettings: mvec " + cat
                        + ":states" + wave, "is not the same size as"
                        " pvec " + names[i]);
      valid = false;
    }
  }

}

//--------------------------------------------------------------------------

// Initialise and check a group of FVec settings.

void SigmaOniaSetup::initSettings(string wave, unsigned int size,
  const vector<string> &names, vector< vector<bool> > &fvecs,
  bool &valid) {

  for (unsigned int i = 0; i < names.size(); ++i) {
    fvecs.push_back(settingsPtr->fvec(names[i]));
    if (fvecs.back().size() != size) {
      infoPtr->errorMsg("Error in SigmaOniaSetup::initSettings: mvec " + cat
                        + ":states" + wave, "is not the same size as"
                        " fvec " + names[i]);
      valid = false;
    }
  }

}

//==========================================================================

// Sigma2gg2QQbar3S11g class.
// Cross section g g -> QQbar[3S1(1)] g (Q = c or b).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2QQbar3S11g::initProc() {

  // Process name.
  nameSave = "g g -> "
    + string((codeSave - codeSave%100)/100 == 4 ? "ccbar" : "bbbar")
    + "(3S1)[3S1(1)] g";

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2gg2QQbar3S11g::sigmaKin() {

  // Calculate kinematics dependence.
  double stH = sH + tH;
  double tuH = tH + uH;
  double usH = uH + sH;
  double sig = (10. * M_PI / 81.) * m3 * ( pow2(sH * tuH)
    + pow2(tH * usH) + pow2(uH * stH) ) / pow2( stH * tuH * usH );

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2QQbar3S11g::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 21);

  // Two orientations of colour flow.
  setColAcol( 1, 2, 2, 3, 0, 0, 1, 3);
  if (rndmPtr->flat() > 0.5) swapColAcol();

}

//==========================================================================

// Sigma2gg2QQbar3S11gm class.
// Cross section g g -> QQbar[3S1(1)] gamma (Q = c or b).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2QQbar3S11gm::initProc() {

  // Process name.
  nameSave = "g g -> "
    + string((codeSave - codeSave%100)/100 == 4 ? "ccbar" : "bbbar")
    + "(3S1)[3S1(1)] gamma";

  // Squared quark charge.
  qEM2 = particleDataPtr->charge((codeSave - codeSave%100)/100);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2gg2QQbar3S11gm::sigmaKin() {

  // Calculate kinematics dependence.
  double stH = sH + tH;
  double tuH = tH + uH;
  double usH = uH + sH;
  double sig = (8. * M_PI / 27.) * m3 * ( pow2(sH * tuH)
    + pow2(tH * usH) + pow2(uH * stH) ) / pow2( stH * tuH * usH );

  // Answer.
  sigma = (M_PI/sH2) * alpEM * qEM2 * pow2(alpS) * oniumME * sig;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2QQbar3S11gm::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 22);

  // Single colour flow orientation.
  setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

}

//==========================================================================

// Sigma2gg2QQbar3PJ1g class.
// Cross section g g -> QQbar[3PJ(1)] g (Q = c or b, J = 0, 1 or 2).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2QQbar3PJ1g::initProc() {

  // Process name.
  if (jSave >= 0 && jSave <= 2)
    nameSave = namePrefix() + " -> " + nameMidfix() + "(3PJ)[3PJ(1)] "
      + namePostfix();
  else
    nameSave = "illegal process";

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2gg2QQbar3PJ1g::sigmaKin() {

  // Useful derived kinematics quantities.
  double pRat  = (sH * uH + uH * tH + tH * sH)/ sH2;
  double qRat  = tH * uH / sH2;
  double rRat  = s3 / sH;
  double pRat2 = pRat * pRat;
  double pRat3 = pRat2 * pRat;
  double pRat4 = pRat3 * pRat;
  double qRat2 = qRat * qRat;
  double qRat3 = qRat2 * qRat;
  double qRat4 = qRat3 * qRat;
  double rRat2 = rRat * rRat;
  double rRat3 = rRat2 * rRat;
  double rRat4 = rRat3 * rRat;

  // Calculate kinematics dependence.
  double sig = 0.;
  if (jSave == 0) {
    sig = (8. * M_PI / (9. * m3 * sH))
      * ( 9. * rRat2 * pRat4 * (rRat4 - 2. * rRat2 * pRat + pRat2)
      - 6. * rRat * pRat3 * qRat * (2. * rRat4 - 5. * rRat2 * pRat
      + pRat2) - pRat2 * qRat2 * (rRat4 + 2. * rRat2 * pRat - pRat2)
      + 2. * rRat * pRat * qRat3 * (rRat2 - pRat) + 6. * rRat2 * qRat4)
      / (qRat * pow4(qRat - rRat * pRat));
  } else if (jSave == 1) {
    sig =  (8. * M_PI / (3.* m3 * sH)) * pRat2
      * (rRat * pRat2 * (rRat2 - 4. * pRat)
      + 2. * qRat * (-rRat4 + 5. * rRat2 * pRat + pRat2)
      - 15. * rRat * qRat2) / pow4(qRat - rRat * pRat);
  } else if (jSave == 2) {
    sig = (8. * M_PI / (9. * m3 * sH))
      * (12. * rRat2 * pRat4 * (rRat4 - 2. * rRat2 * pRat + pRat2)
      - 3. * rRat * pRat3 * qRat * (8. * rRat4 - rRat2 * pRat + 4. * pRat2)
      + 2. * pRat2 * qRat2 * (-7. * rRat4 + 43. * rRat2 * pRat + pRat2)
      + rRat * pRat * qRat3 * (16. * rRat2 - 61. * pRat)
      + 12. * rRat2 * qRat4) / (qRat * pow4(qRat-rRat * pRat));
  }

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2QQbar3PJ1g::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 21);

  // Two orientations of colour flow.
  setColAcol( 1, 2, 2, 3, 0, 0, 1, 3);
  if (rndmPtr->flat() > 0.5) swapColAcol();

}

//==========================================================================

// Sigma2qg2QQbar3PJ1q class.
// Cross section q g -> QQbar[3PJ(1)] q (Q = c or b, J = 0, 1 or 2).

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2qg2QQbar3PJ1q::sigmaKin() {

  // Calculate kinematics dependence.
  double usH = uH + sH;
  double sig = 0.;
  if (jSave == 0) {
    sig = - (16. * M_PI / 81.) * pow2(tH - 3. * s3) * (sH2 + uH2)
      / (m3 * tH * pow4(usH));
  } else if (jSave == 1) {
    sig = - (32. * M_PI / 27.) * (4. * s3 * sH * uH + tH * (sH2 + uH2))
      / (m3 * pow4(usH));
  } else if (jSave == 2) {
    sig = - (32. *M_PI / 81.) * ( (6. * s3*s3 + tH2) * pow2(usH)
      - 2. * sH * uH * (tH2 + 6. * s3 * usH)) / (m3 * tH * pow4(usH));
  }

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2QQbar3PJ1q::setIdColAcol() {

  // Flavours are trivial.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, idHad, idq);

  // tH defined between q_in and q_out: must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21);

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2QQbar3PJ1g class.
// Cross section q qbar -> QQbar[3PJ(1)] g (Q = c or b, J = 0, 1 or 2).

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2qqbar2QQbar3PJ1g::sigmaKin() {

  // Calculate kinematics dependence.
  double tuH = tH + uH;
  double sig = 0.;
  if (jSave == 0) {
    sig =(128. * M_PI / 243.) * pow2(sH - 3. * s3) * (tH2 + uH2)
      / (m3 * sH * pow4(tuH));
  } else if (jSave == 1) {
    sig = (256. * M_PI / 81.) * (4. * s3 * tH * uH + sH * (tH2 + uH2))
      / (m3 * pow4(tuH));
  } else if (jSave == 2) {
    sig = (256. * M_PI / 243.) * ( (6. * s3*s3 + sH2) * pow2(tuH)
      - 2. * tH * uH * (sH2 + 6. * s3 * tuH) )/ (m3 * sH * pow4(tuH));
  }

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2QQbar3PJ1g::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 21);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2gg2QQbar3DJ1g class.
// Cross section g g -> QQbar[3DJ(1)] g (Q = c or b).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2QQbar3DJ1g::initProc() {

  // Process name.
  if (jSave >= 1 && jSave <= 3)
    nameSave = namePrefix() + " -> " + nameMidfix() + "(3DJ)[3DJ(1)] "
      + namePostfix();
  else
    nameSave = "illegal process";

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2gg2QQbar3DJ1g::sigmaKin() {

  // Calculate kinematics dependence.
  double m2V[12], sHV[12], mpsV[8], mmsV[6], mmtV[6], sptV[6];
  m2V[0]  = 1;
  sHV[0]  = 1;
  mmtV[0] = 1;
  mpsV[0] = 1;
  mmsV[0] = 1;
  sptV[0] = 1;
  for (int i = 1; i < 12; ++i) {
    m2V[i] = m2V[i - 1] * s3;
    sHV[i] = sHV[i - 1] * sH;
    if (i < 8) {
      mpsV[i] = mpsV[i - 1] * (s3 + sH);
      if (i < 6) {
        mmsV[i] = mmsV[i - 1] * (s3 - sH);
        mmtV[i] = mmtV[i - 1] * (s3 - tH);
        sptV[i] = sptV[i - 1] * (sH + tH);
      }
    }
  }
  double fac = (pow3(alpS)*pow2(M_PI));
  double sig = 0;
  if (jSave == 1) {
    fac *= 16. / 81.;
    sig  = -25/(sqrt(m2V[1])*mmsV[5]) + (49*sqrt(m2V[3]))/(mmsV[5]*sHV[2])
      + (48*sqrt(m2V[3])*sHV[2]*(m2V[2] + sHV[2]))/(mmsV[3]*mmtV[5]*mpsV[3])
      - (67*sqrt(m2V[1]))/(mmsV[5]*sHV[1]) - (5*sHV[1])/(sqrt(m2V[3])*mmsV[5])
      + (4*sqrt(m2V[1])*(m2V[6] + 97*m2V[4]*sHV[2] - 48*m2V[3]*sHV[3]
      + 105*m2V[2]*sHV[4] + 33*sHV[6] -
      24*m2V[5]*sHV[1]))/(mmsV[4]*mmtV[4]*mpsV[4]) - (4*(m2V[9] +
      197*m2V[7]*sHV[2] - 50*m2V[6]*sHV[3] + 509*m2V[5]*sHV[4] -
      416*m2V[4]*sHV[5] + 237*m2V[3]*sHV[6] - 400*m2V[2]*sHV[7] -
      10*sHV[9] -
      164*m2V[8]*sHV[1]))/(sqrt(m2V[1])*mmsV[5]*mmtV[3]*mpsV[5]*sHV[1]) +
      (224*m2V[10] + 1825*m2V[8]*sHV[2] - 3980*m2V[7]*sHV[3] +
      3996*m2V[6]*sHV[4] - 4766*m2V[5]*sHV[5] + 10022*m2V[4]*sHV[6] -
      5212*m2V[3]*sHV[7] + 6124*m2V[2]*sHV[8] - 869*m2V[1]*sHV[9] +
      145*sHV[10] -
      597*m2V[9]*sHV[1])/(sqrt(m2V[1])*mmsV[5]*mmtV[1]*mpsV[7]*sHV[2]) +
      (102*m2V[11] + 331*m2V[9]*sHV[2] - 2021*m2V[8]*sHV[3] +
      3616*m2V[7]*sHV[4] - 968*m2V[6]*sHV[5] + 3386*m2V[5]*sHV[6] -
      6150*m2V[4]*sHV[7] + 666*m2V[3]*sHV[8] - 1134*m2V[2]*sHV[9] -
      5*m2V[1]*sHV[10] - 5*sHV[11] -
      506*m2V[10]*sHV[1])/(sqrt(m2V[3])*mmsV[5]*mmtV[2]*mpsV[6]*sHV[2]) +
      (48*sqrt(m2V[3])*sHV[2]*(m2V[2] + sHV[2]))/(mmsV[3]*mpsV[3]*sptV[5])
      + (4*sqrt(m2V[1])*(m2V[6] + 97*m2V[4]*sHV[2] - 48*m2V[3]*sHV[3] +
      105*m2V[2]*sHV[4] + 33*sHV[6] -
      24*m2V[5]*sHV[1]))/(mmsV[4]*mpsV[4]*sptV[4]) - (4*(m2V[9] +
      197*m2V[7]*sHV[2] - 50*m2V[6]*sHV[3] + 509*m2V[5]*sHV[4] -
      416*m2V[4]*sHV[5] + 237*m2V[3]*sHV[6] - 400*m2V[2]*sHV[7] -
      10*sHV[9] -
      164*m2V[8]*sHV[1]))/(sqrt(m2V[1])*mmsV[5]*mpsV[5]*sHV[1]*sptV[3]) +
      (102*m2V[11] + 331*m2V[9]*sHV[2] - 2021*m2V[8]*sHV[3] +
      3616*m2V[7]*sHV[4] - 968*m2V[6]*sHV[5] + 3386*m2V[5]*sHV[6] -
      6150*m2V[4]*sHV[7] + 666*m2V[3]*sHV[8] - 1134*m2V[2]*sHV[9] -
      5*m2V[1]*sHV[10] - 5*sHV[11] -
      506*m2V[10]*sHV[1])/(sqrt(m2V[3])*mmsV[5]*mpsV[6]*sHV[2]*sptV[2]) +
      (224*m2V[10] + 1825*m2V[8]*sHV[2] - 3980*m2V[7]*sHV[3] +
      3996*m2V[6]*sHV[4] - 4766*m2V[5]*sHV[5] + 10022*m2V[4]*sHV[6] -
      5212*m2V[3]*sHV[7] + 6124*m2V[2]*sHV[8] - 869*m2V[1]*sHV[9] +
      145*sHV[10] -
      597*m2V[9]*sHV[1])/(sqrt(m2V[1])*mmsV[5]*mpsV[7]*sHV[2]*sptV[1]);
  } else if (jSave == 2) {
    fac *= 32. / 27.;
    sig  = 16/(sqrt(m2V[1])*mmsV[5]) +
      (2*sqrt(m2V[3]))/(mmsV[5]*sHV[2]) - (8*sqrt(m2V[3])*sHV[2]*(m2V[2] +
      sHV[2]))/(mmsV[3]*mmtV[5]*mpsV[3]) +
      (6*sqrt(m2V[1]))/(mmsV[5]*sHV[1]) -
      (16*sHV[1])/(sqrt(m2V[3])*mmsV[5]) - (2*sqrt(m2V[1])*(3*m2V[6] -
      25*m2V[4]*sHV[2] - 16*m2V[3]*sHV[3] - 33*m2V[2]*sHV[4] - 5*sHV[6] -
      8*m2V[5]*sHV[1]))/(mmsV[4]*mmtV[4]*mpsV[4]) + (2*(3*m2V[9] -
      41*m2V[7]*sHV[2] - 37*m2V[6]*sHV[3] - 149*m2V[5]*sHV[4] +
      55*m2V[4]*sHV[5] - 53*m2V[3]*sHV[6] + 167*m2V[2]*sHV[7] + 16*sHV[9]
      + 7*m2V[8]*sHV[1]))/(sqrt(m2V[1])*mmsV[5]*mmtV[3]*mpsV[5]*sHV[1]) +
      (2*(m2V[10] + 34*m2V[8]*sHV[2] - 198*m2V[7]*sHV[3] -
      140*m2V[6]*sHV[4] - 746*m2V[5]*sHV[5] + 226*m2V[4]*sHV[6] -
      486*m2V[3]*sHV[7] + 679*m2V[2]*sHV[8] - 50*m2V[1]*sHV[9] +
      112*sHV[10] -
      8*m2V[9]*sHV[1]))/(sqrt(m2V[1])*mmsV[5]*mmtV[1]*mpsV[7]*sHV[2]) +
      (m2V[11] + 19*m2V[9]*sHV[2] - m2V[8]*sHV[3] + 597*m2V[7]*sHV[4] +
      321*m2V[6]*sHV[5] + 797*m2V[5]*sHV[6] - 791*m2V[4]*sHV[7] +
      26*m2V[3]*sHV[8] - 468*m2V[2]*sHV[9] - 16*m2V[1]*sHV[10] -
      16*sHV[11] -
      21*m2V[10]*sHV[1])/(sqrt(m2V[3])*mmsV[5]*mmtV[2]*mpsV[6]*sHV[2]) -
      (8*sqrt(m2V[3])*sHV[2]*(m2V[2] + sHV[2]))/(mmsV[3]*mpsV[3]*sptV[5])
      - (2*sqrt(m2V[1])*(3*m2V[6] - 25*m2V[4]*sHV[2] - 16*m2V[3]*sHV[3] -
      33*m2V[2]*sHV[4] - 5*sHV[6] -
      8*m2V[5]*sHV[1]))/(mmsV[4]*mpsV[4]*sptV[4]) + (2*(3*m2V[9] -
      41*m2V[7]*sHV[2] - 37*m2V[6]*sHV[3] - 149*m2V[5]*sHV[4] +
      55*m2V[4]*sHV[5] - 53*m2V[3]*sHV[6] + 167*m2V[2]*sHV[7] + 16*sHV[9]
      + 7*m2V[8]*sHV[1]))/(sqrt(m2V[1])*mmsV[5]*mpsV[5]*sHV[1]*sptV[3]) +
      (m2V[11] + 19*m2V[9]*sHV[2] - m2V[8]*sHV[3] + 597*m2V[7]*sHV[4] +
      321*m2V[6]*sHV[5] + 797*m2V[5]*sHV[6] - 791*m2V[4]*sHV[7] +
      26*m2V[3]*sHV[8] - 468*m2V[2]*sHV[9] - 16*m2V[1]*sHV[10] -
      16*sHV[11] -
      21*m2V[10]*sHV[1])/(sqrt(m2V[3])*mmsV[5]*mpsV[6]*sHV[2]*sptV[2]) +
      (2*(m2V[10] + 34*m2V[8]*sHV[2] - 198*m2V[7]*sHV[3] -
      140*m2V[6]*sHV[4] - 746*m2V[5]*sHV[5] + 226*m2V[4]*sHV[6] -
      486*m2V[3]*sHV[7] + 679*m2V[2]*sHV[8] - 50*m2V[1]*sHV[9] +
      112*sHV[10] -
      8*m2V[9]*sHV[1]))/(sqrt(m2V[1])*mmsV[5]*mpsV[7]*sHV[2]*sptV[1]);
  } else if (jSave == 3) {
    fac *= 256. / 189.;
    sig  = 5/(sqrt(m2V[1])*mmsV[5]) + sqrt(m2V[3])/(mmsV[5]*sHV[2]) +
      (2*sqrt(m2V[3])*sHV[2]*(m2V[2] + sHV[2]))/(mmsV[3]*mmtV[5]*mpsV[3])
      - (3*sqrt(m2V[1]))/(mmsV[5]*sHV[1]) -
      (5*sHV[1])/(sqrt(m2V[3])*mmsV[5]) + (sqrt(m2V[1])*(6*m2V[6] +
      67*m2V[4]*sHV[2] - 8*m2V[3]*sHV[3] + 45*m2V[2]*sHV[4] + 8*sHV[6] -
      4*m2V[5]*sHV[1]))/(mmsV[4]*mmtV[4]*mpsV[4]) + (-6*m2V[9] -
      152*m2V[7]*sHV[2] + 80*m2V[6]*sHV[3] - 269*m2V[5]*sHV[4] +
      211*m2V[4]*sHV[5] - 77*m2V[3]*sHV[6] + 155*m2V[2]*sHV[7] + 10*sHV[9]
      + 64*m2V[8]*sHV[1])/(sqrt(m2V[1])*mmsV[5]*mmtV[3]*mpsV[5]*sHV[1]) +
      (16*m2V[10] + 295*m2V[8]*sHV[2] - 555*m2V[7]*sHV[3] +
      769*m2V[6]*sHV[4] - 1079*m2V[5]*sHV[5] + 913*m2V[4]*sHV[6] -
      603*m2V[3]*sHV[7] + 601*m2V[2]*sHV[8] - 56*m2V[1]*sHV[9] +
      70*sHV[10] -
      83*m2V[9]*sHV[1])/(sqrt(m2V[1])*mmsV[5]*mmtV[1]*mpsV[7]*sHV[2]) +
      (8*m2V[11] + 104*m2V[9]*sHV[2] - 284*m2V[8]*sHV[3] +
      549*m2V[7]*sHV[4] - 282*m2V[6]*sHV[5] + 514*m2V[5]*sHV[6] -
      520*m2V[4]*sHV[7] + 34*m2V[3]*sHV[8] - 171*m2V[2]*sHV[9] -
      5*m2V[1]*sHV[10] - 5*sHV[11] -
      54*m2V[10]*sHV[1])/(sqrt(m2V[3])*mmsV[5]*mmtV[2]*mpsV[6]*sHV[2]) +
      (2*sqrt(m2V[3])*sHV[2]*(m2V[2] + sHV[2]))/(mmsV[3]*mpsV[3]*sptV[5])
      + (sqrt(m2V[1])*(6*m2V[6] + 67*m2V[4]*sHV[2] - 8*m2V[3]*sHV[3] +
      45*m2V[2]*sHV[4] + 8*sHV[6] -
      4*m2V[5]*sHV[1]))/(mmsV[4]*mpsV[4]*sptV[4]) + (-6*m2V[9] -
      152*m2V[7]*sHV[2] + 80*m2V[6]*sHV[3] - 269*m2V[5]*sHV[4] +
      211*m2V[4]*sHV[5] - 77*m2V[3]*sHV[6] + 155*m2V[2]*sHV[7] + 10*sHV[9]
      + 64*m2V[8]*sHV[1])/(sqrt(m2V[1])*mmsV[5]*mpsV[5]*sHV[1]*sptV[3]) +
      (8*m2V[11] + 104*m2V[9]*sHV[2] - 284*m2V[8]*sHV[3] +
      549*m2V[7]*sHV[4] - 282*m2V[6]*sHV[5] + 514*m2V[5]*sHV[6] -
      520*m2V[4]*sHV[7] + 34*m2V[3]*sHV[8] - 171*m2V[2]*sHV[9] -
      5*m2V[1]*sHV[10] - 5*sHV[11] -
      54*m2V[10]*sHV[1])/(sqrt(m2V[3])*mmsV[5]*mpsV[6]*sHV[2]*sptV[2]) +
      (16*m2V[10] + 295*m2V[8]*sHV[2] - 555*m2V[7]*sHV[3] +
      769*m2V[6]*sHV[4] - 1079*m2V[5]*sHV[5] + 913*m2V[4]*sHV[6] -
      603*m2V[3]*sHV[7] + 601*m2V[2]*sHV[8] - 56*m2V[1]*sHV[9] +
      70*sHV[10] -
      83*m2V[9]*sHV[1])/(sqrt(m2V[1])*mmsV[5]*mpsV[7]*sHV[2]*sptV[1]);
  }

  // Answer.
  sigma = ((2.*jSave + 1.) / 3.) * oniumME * fac * sig;

}

//==========================================================================

// Sigma2gg2QQbarX8g class.
// Cross section g g -> QQbar[X(8)] g (Q = c or b, X = 3S1, 1S0 or 3PJ).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2QQbarX8g::initProc() {

  // Return for illegal process.
  if (stateSave < 0 || stateSave > 2) {
    idHad = 0;
    nameSave = "illegal process";
    return;
  }

  // Determine quark composition and quantum numbers.
  int mod1(10), mod2(1);
  vector<int> digits;
  while (digits.size() < 7) {
    digits.push_back((idHad%mod1 - idHad%mod2) / mod2);
    mod1 *= 10;
    mod2 *= 10;
  }
  int s, l, j((digits[0] - 1)/2);
  if (j != 0) {
    if      (digits[4] == 0) {l = j - 1; s = 1;}
    else if (digits[4] == 1) {l = j;     s = 0;}
    else if (digits[4] == 2) {l = j;     s = 1;}
    else                     {l = j + 1; s = 1;}
  } else {
    if      (digits[4] == 0) {l = 0;  s = 0;}
    else                     {l = 1;  s = 1;}
  }

  // Set the process name.
  stringstream sName, jName;
  string lName, stateName;
  sName << 2*s + 1;
  if (l == 0) jName << j;
  else jName << "J";
  if (l == 0) lName = "S";
  else if (l == 1) lName = "P";
  else if (l == 2) lName = "D";
  if (stateSave == 0) stateName = "[3S1(8)]";
  else if (stateSave == 1) stateName = "[1S0(8)]";
  else if (stateSave == 2) stateName = "[3PJ(8)]";
  nameSave = namePrefix() + " -> " + (digits[1] == 4 ? "ccbar" : "bbbar")
    + "(" + sName.str() + lName + jName.str() + ")" + stateName
    + " " + namePostfix();

  // Ensure the dummy particle for the colour-octet state is valid.
  int idOct = 9900000 + digits[1]*10000 + stateSave*1000 + digits[5]*100
    + digits[4]*10 + digits[0];
  double m0     = particleDataPtr->m0(idHad) + abs(mSplit);
  double mWidth = 0.0;
  if (!particleDataPtr->isParticle(idOct)) {
    string nameOct    = particleDataPtr->name(idHad) + stateName;
    int    spinType   = stateSave == 1 ? 1 : 3;
    int    chargeType = particleDataPtr->chargeType(idHad);
    int    colType    = 2;
    particleDataPtr->addParticle(idOct, nameOct, spinType, chargeType, colType,
                                 m0, mWidth, m0, m0);
    ParticleDataEntry* entry = particleDataPtr->particleDataEntryPtr(idOct);
    if (entry) entry->addChannel(1, 1.0, 0, idHad, 21);
  } else if (mSplit > 0 && abs(particleDataPtr->m0(idOct) - m0) > 1E-5) {
    particleDataPtr->m0(idOct, m0);
    particleDataPtr->mWidth(idOct, mWidth);
    particleDataPtr->mMin(idOct, m0);
    particleDataPtr->mMax(idOct, m0);
  } else if (particleDataPtr->m0(idOct) <= particleDataPtr->m0(idHad)) {
    infoPtr->errorMsg("Warning in Sigma2gg2QQbarX8g::initProc: mass of "
                      "intermediate colour-octet state"
                      "increased to be greater than the physical state");
    particleDataPtr->m0(idOct, m0);
    particleDataPtr->mWidth(idOct, mWidth);
    particleDataPtr->mMin(idOct, m0);
    particleDataPtr->mMax(idOct, m0);
  }
  idHad = idOct;

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2gg2QQbarX8g::sigmaKin() {

  // Calculate kinematics dependence.
  double stH = sH + tH;
  double tuH = tH + uH;
  double usH = uH + sH;
  double sig = 0.;
  if (stateSave == 0) {
    sig = (M_PI / 72.) * m3 * ( 27. * (pow2(stH) + pow2(tuH)
      + pow2(usH)) / (s3*s3) - 16. ) * ( pow2(sH * tuH)
      + pow2(tH * usH) + pow2(uH * stH) ) / pow2( stH * tuH * usH );
  } else if (stateSave == 1) {
    sig = (5. * M_PI / 16.) * m3 * ( pow2(uH / (tuH * usH))
      + pow2(sH / (stH * usH)) + pow2(tH / (stH * tuH)) ) * ( 12.
      + (pow4(stH) + pow4(tuH) + pow4(usH)) / (s3 * sH * tH * uH) );
  } else if (stateSave == 2) {
    double sH3 = sH2 * sH;
    double sH4 = sH3 * sH;
    double sH5 = sH4 * sH;
    double sH6 = sH5 * sH;
    double sH7 = sH6 * sH;
    double sH8 = sH7 * sH;
    double tH3 = tH2 * tH;
    double tH4 = tH3 * tH;
    double tH5 = tH4 * tH;
    double tH6 = tH5 * tH;
    double tH7 = tH6 * tH;
    double tH8 = tH7 * tH;
    double ssttH = sH * sH + sH * tH + tH * tH;
    sig = 5. * M_PI * (3. * sH * tH * stH * pow4(ssttH)
      - s3 * pow2(ssttH) * (7. * sH6 + 36. * sH5 * tH + 45. * sH4 * tH2
        + 28. * sH3 * tH3 + 45. * sH2 * tH4 + 36. * sH * tH5 + 7. * tH6)
      + pow2(s3) * stH * (35. *sH8 + 169. * sH7 * tH + 299. * sH6 * tH2
        + 401. * sH5 * tH3 + 418. * sH4 * tH4 + 401. * sH3 * tH5
        + 299. * sH2 * tH6 + 169. * sH * tH7 + 35. * tH8)
      - pow3(s3) * (84. *sH8+432. *sH7*tH+905. *sH6*tH2
        + 1287. * sH5 * tH3 + 1436. * sH4 * tH4 +1287. * sH3 * tH5
        + 905. * sH2 * tH6 + 432. * sH * tH7 + 84. * tH8)
      + pow4(s3) * stH * (126. * sH6 + 451. * sH5 * tH +677. * sH4 * tH2
        + 836. * sH3 * tH3 + 677. * sH2 * tH4 + 451. * sH * tH5
        + 126. * tH6)
      - pow5(s3) * 3. * (42. * sH6 + 171. * sH5 * tH + 304. * sH4 * tH2
        + 362. * sH3 * tH3 + 304. * sH2 * tH4 + 171. * sH * tH5
        + 42. * tH6)
      + pow3(s3 * s3) * 2. * stH * (42. * sH4 + 106. * sH3 * tH
        + 119. * sH2 * tH2 + 106. * sH * tH3 + 42. * tH4)
      - pow4(s3) * pow3(s3) * (35. * sH4 + 99. * sH3 * tH
        + 120. * sH2 * tH2 + 99.  * sH * tH3 + 35.  * tH4)
      + pow4(s3 * s3) * 7. * stH * ssttH)
      / (sH * tH * uH * s3 * m3 * pow3(stH * tuH * usH));
  }

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2QQbarX8g::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 21);

  // Split total contribution into different colour flows just like in
  // g g -> g g (with kinematics recalculated for massless partons).
  double sHr    = - (tH + uH);
  double sH2r   = sHr * sHr;
  double sigTS  = tH2/sH2r + 2.*tH/sHr + 3. + 2.*sHr/tH + sH2r/tH2;
  double sigUS  = uH2/sH2r + 2.*uH/sHr + 3. + 2.*sHr/uH + sH2r/uH2;
  double sigTU  = tH2/uH2 + 2.*tH/uH + 3. + 2.*uH/tH + uH2/tH2;
  double sigSum = sigTS + sigUS + sigTU;

  // Three colour flow topologies, each with two orientations.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 4, 4, 3);
  else if (sigRand < sigTS + sigUS)
                       setColAcol( 1, 2, 3, 1, 3, 4, 4, 2);
  else                 setColAcol( 1, 2, 3, 4, 1, 4, 3, 2);
  if (rndmPtr->flat() > 0.5) swapColAcol();

}

//==========================================================================

// Sigma2qg2QQbarX8q class.
// Cross section q g -> QQbar[X(8)] q (Q = c or b, X = 3S1, 1S0 or 3PJ).

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2qg2QQbarX8q::sigmaKin() {

  // Calculate kinematics dependence.
  double stH  = sH + tH;
  double tuH  = tH + uH;
  double usH  = uH + sH;
  double stH2 = stH * stH;
  double tuH2 = tuH * tuH;
  double usH2 = usH * usH;
  double sig  = 0.;
  if (stateSave == 0) {
    sig = - (M_PI / 27.)* (4. * (sH2 + uH2) - sH * uH) * (stH2 +tuH2)
      / (s3 * m3 * sH * uH * usH2);
  } else if (stateSave == 1) {
    sig = - (5. * M_PI / 18.) * (sH2 + uH2) / (m3 * tH * usH2);
  } else if (stateSave == 2) {
    sig = - (10. * M_PI / 9.) * ( (7. * usH + 8. * tH) * (sH2 + uH2)
      + 4. * tH * (2. * pow2(s3) - stH2 - tuH2) )
      / (s3 * m3 * tH * usH2 * usH);
  }

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qg2QQbarX8q::setIdColAcol() {

  // Flavours are trivial.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, idHad, idq);

  // tH defined between q_in and q_out: must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21);

  // Split total contribution into different colour flows just like in
  // q g -> q g (with kinematics recalculated for massless partons).
  double sHr    = - (tH + uH);
  double sH2r   = sHr * sHr;
  double sigTS  = uH2/tH2 - (4./9.) * uH/sHr;
  double sigTU  = sH2r/tH2 - (4./9.) * sHr/uH;
  double sigSum = sigTS + sigTU;

  // Two colour flow topologies. Swap if first is gluon, or when antiquark.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 2, 1, 2, 3, 3, 0);
  else                 setColAcol( 1, 0, 2, 3, 1, 3, 2, 0);
  if (id1 == 21) swapCol12();
  if (idq < 0) swapColAcol();

}

//==========================================================================

// Sigma2qqbar2QQbarX8g class.
// Cross section q qbar -> QQbar[X(8)] g (Q = c or b, X = 3S1, 1S0 or 3PJ).

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2qqbar2QQbarX8g::sigmaKin() {

  // Calculate kinematics dependence.
  double stH  = sH + tH;
  double tuH  = tH + uH;
  double usH  = uH + sH;
  double stH2 = stH * stH;
  double tuH2 = tuH * tuH;
  double usH2 = usH * usH;
  double sig  = 0.;
  if (stateSave == 0) {
    sig = (8. * M_PI / 81.) * (4. * (tH2 + uH2) - tH * uH)
      * (stH2 + usH2) / (s3 * m3 * tH * uH * tuH2);
  } else if (stateSave == 1) {
    sig = (20. * M_PI / 27.) * (tH2 + uH2) / (m3 * sH * tuH2);
  } else if (stateSave == 2) {
    sig = (80. * M_PI / 27.) * ( (7. * tuH + 8. * sH) * (tH2 + uH2)
      + 4. * sH * (2. * pow2(s3) - stH2 -usH2) )
      / (s3 * m3 * sH * tuH2 * tuH);
  }

  // Answer.
  sigma = (M_PI/sH2) * pow3(alpS) * oniumME * sig;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2QQbarX8g::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad, 21);

  // Split total contribution into different colour flows just like in
  // q qbar -> g g (with kinematics recalculated for massless partons).
  double sHr    = - (tH + uH);
  double sH2r   = sHr * sHr;
  double sigTS  = (4. / 9.) * uH / tH - uH2 / sH2r;
  double sigUS  = (4. / 9.) * tH / uH - tH2 / sH2r;
  double sigSum = sigTS + sigUS;

  // Two colour flow topologies. Swap if first is antiquark.
  double sigRand = sigSum * rndmPtr->flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 0, 2, 1, 3, 3, 2);
  else                 setColAcol( 1, 0, 0, 2, 3, 2, 1, 3);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

// Sigma2gg2QQbar3S11QQbar3S11 class.
// Cross section g g -> QQbar[3S1(1)] QQbar[3S1(1)] (Q = c or b).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2gg2QQbar3S11QQbar3S11::initProc() {

  // Process name.
  int flavor((codeSave - codeSave%100)/100);
  nameSave = string(flavor == 4 ? "ccbar" : "bbbar");
  nameSave = "g g -> double " + nameSave + "(3S1)[3S1(1)]";

  // Constant mass squared vector.
  m2V.push_back(1.0); m2V.push_back(pow2(2. * particleDataPtr->m0(flavor)));
  for (int iSqr = 2; iSqr < 14; ++iSqr) m2V.push_back(m2V[iSqr - 1] * m2V[1]);

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2gg2QQbar3S11QQbar3S11::sigmaKin() {

  // Values of sH, uH, and tH with exponents above 6.
  double tH7(pow6(tH) * tH),tH8(tH7 * tH), tH9(tH8 * tH), tH10(tH9 * tH),
    uH7(pow6(uH) * uH),uH8(uH7 * uH), uH9(uH8 * uH), uH10(uH9 * uH),
    sH8(pow6(sH) * sH * sH);

  // The answer.
  sigma = (64*pow4(alpS)*oniumME1*oniumME2*pow3(M_PI)*(2680*m2V[12]
    - 14984*m2V[11]*(tH + uH) - 16*m2V[9]*(tH + uH)*(1989*pow2(tH)
    + 10672*tH*uH + 1989*pow2(uH)) + m2V[10]*(31406*pow2(tH) + 89948*tH*uH
    + 31406*pow2(uH)) + 2*pow4(tH)*pow4(uH)*(349*pow4(tH) - 908*pow3(tH)*uH
    + 1374*pow2(tH)*pow2(uH) - 908*tH*pow3(uH) + 349*pow4(uH)) - 4*m2V[7]*(tH
    + uH)*(1793*pow4(tH) + 36547*pow3(tH)*uH + 97572*pow2(tH)*pow2(uH)
    + 36547*tH*pow3(uH) + 1793*pow4(uH)) + 4*m2V[8]*(4417*pow4(tH)
    + 57140*pow3(tH)*uH + 117714*pow2(tH)*pow2(uH) + 57140*tH*pow3(uH)
    + 4417*pow4(uH)) + 4*m2V[1]*pow2(tH)*pow2(uH)*(tH + uH)*(9*pow6(tH)
    - 595*pow5(tH)*uH + 558*pow4(tH)*pow2(uH) - 952*pow3(tH)*pow3(uH)
    + 558*pow2(tH)*pow4(uH) - 595*tH*pow5(uH) + 9*pow6(uH)) - 2*m2V[5]*(tH
    + uH)*(397*pow6(tH) + 14994*pow5(tH)*uH + 76233*pow4(tH)*pow2(uH)
    + 91360*pow3(tH)*pow3(uH) + 76233*pow2(tH)*pow4(uH) + 14994*tH*pow5(uH)
    + 397*pow6(uH)) + m2V[6]*(2956*pow6(tH) + 76406*pow5(tH)*uH
    + 361624*pow4(tH)*pow2(uH) + 571900*pow3(tH)*pow3(uH)
    + 361624*pow2(tH)*pow4(uH) + 76406*tH*pow5(uH) + 2956*pow6(uH))
    + 2*m2V[3]*(tH + uH)*(10*tH8 - 421*tH7*uH - 8530*pow6(tH)*pow2(uH)
    - 20533*pow5(tH)*pow3(uH) + 2880*pow4(tH)*pow4(uH)
    - 20533*pow3(tH)*pow5(uH) - 8530*pow2(tH)*pow6(uH) - 421*tH*uH7 + 10*uH8)
    + m2V[4]*(47*tH8 + 7642*tH7*uH + 73146*pow6(tH)*pow2(uH)
    + 150334*pow5(tH)*pow3(uH) + 132502*pow4(tH)*pow4(uH)
    + 150334*pow3(tH)*pow5(uH) + 73146*pow2(tH)*pow6(uH) + 7642*tH*uH7
    + 47*uH8) + m2V[2]*(tH10 - 66*tH9*uH + 2469*tH8*pow2(uH)
    + 12874*tH7*pow3(uH) + 11928*pow6(tH)*pow4(uH) + 1164*pow5(tH)*pow5(uH)
    + 11928*pow4(tH)*pow6(uH) + 12874*pow3(tH)*uH7 + 2469*pow2(tH)*uH8
    - 66*tH*uH9 + uH10)))/(6561.*m2V[1]*sH8*pow4(m2V[1] - tH)*pow4(m2V[1]
    - uH));
  if (idHad1 != idHad2) sigma *= 2.;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2gg2QQbar3S11QQbar3S11::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad1, idHad2);

  // One orientation of colour flow.
  setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

}


//==========================================================================

// Sigma2qqbar2QQbar3S11QQbar3S11 class.
// Cross section q qbar -> QQbar[3S1(1)] QQbar[3S1(1)] (Q = c or b).

//--------------------------------------------------------------------------

// Initialize process.

void Sigma2qqbar2QQbar3S11QQbar3S11::initProc() {

  // Process name.
  int flavor((codeSave - codeSave%100)/100);
  nameSave = string(flavor == 4 ? "ccbar" : "bbbar");
  nameSave = "q qbar -> double " + nameSave + "(3S1)[3S1(1)]";

  // Constant mass squared.
  m2 = pow2(2. * particleDataPtr->m0(flavor));

}

//--------------------------------------------------------------------------

// Evaluate d(sigmaHat)/d(tHat); no explicit flavour dependence.

void Sigma2qqbar2QQbar3S11QQbar3S11::sigmaKin() {

  // The answer.
  sigma = (16384*pow4(alpS)*oniumME1*oniumME2*pow3(M_PI)*(6*pow4(sH)
    - 5*pow2(sH)*pow2(tH - uH) - 3*pow4(tH - uH) + 4*pow3(sH)*(tH + uH)
    - 6*sH*pow2(tH - uH)*(tH + uH)))/(19683.*m2*pow6(sH)*pow2(sH));
  if (idHad1 != idHad2) sigma *= 2.;

}

//--------------------------------------------------------------------------

// Select identity, colour and anticolour.

void Sigma2qqbar2QQbar3S11QQbar3S11::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idHad1, idHad2);

  // Two orientations of colour flow.
  setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//==========================================================================

} // end namespace Pythia8
