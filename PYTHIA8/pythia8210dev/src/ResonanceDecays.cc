// ResonanceDecays.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for
// the ResonanceDecays class.

#include "Pythia8/ResonanceDecays.h"

namespace Pythia8 {

//==========================================================================

// The ResonanceDecays class.
// Do all resonance decays sequentially.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of tries to pick a decay channel.
const int    ResonanceDecays::NTRYCHANNEL = 10;

// Number of tries to pick a set of daughter masses.
const int    ResonanceDecays::NTRYMASSES  = 10000;

// Mass above threshold for allowed decays.
const double ResonanceDecays::MSAFETY     = 0.1;

// When constrainted kinematics cut high-mass tail of Breit-Wigner.
const double ResonanceDecays::WIDTHCUT    = 5.;

// Small number (relative to 1) to protect against roundoff errors.
const double ResonanceDecays::TINY        = 1e-10;

// Forbid small Breit-Wigner mass range, as mapped onto atan range.
const double ResonanceDecays::TINYBWRANGE = 1e-8;

// These numbers are hardwired empirical parameters,
// intended to speed up the M-generator.
const double ResonanceDecays::WTCORRECTION[11] = { 1., 1., 1.,
  2., 5., 15., 60., 250., 1250., 7000., 50000. };

//--------------------------------------------------------------------------

bool ResonanceDecays::next( Event& process, int iDecNow) {

  // Loop over all entries to find resonances that should decay.
  // (Except for iDecNow > 0, where only it will be handled.)
  int iDec = iDecNow;
  do {
    Particle& decayer = process[iDec];
    if (decayer.isFinal() && decayer.canDecay() && decayer.mayDecay()
    && decayer.isResonance() ) {

      // Fill the decaying particle in slot 0 of arrays.
      id0    = decayer.id();
      m0     = decayer.m();
      idProd.resize(0);
      mProd.resize(0);
      idProd.push_back( id0 );
      mProd.push_back( m0 );

      // Mother flavour - relevant for gamma*/Z0 mixing. (Not always??)
      int idIn = process[decayer.mother1()].id();

      // Prepare decay selection.
      if (!decayer.particleDataEntry().preparePick(id0, m0, idIn)) {
        ostringstream osWarn;
        osWarn << "for id = " << id0;
        infoPtr->errorMsg("Error in ResonanceDecays::next:"
          " no open decay channel", osWarn.str());
        return false;
      }

      // Pick a decay channel; allow up to ten tries.
      bool foundChannel = false;
      for (int iTryChannel = 0; iTryChannel < NTRYCHANNEL; ++iTryChannel) {

        // Pick decay channel. Find multiplicity.
        DecayChannel& channel = decayer.particleDataEntry().pickChannel();
        mult = channel.multiplicity();

        // Read out flavours.
        idProd.resize(1);
        int idNow;
        for (int i = 1; i <= mult; ++i) {
          idNow = channel.product(i - 1);
          if (id0 < 0 && particleDataPtr->hasAnti(idNow)) idNow = -idNow;
          idProd.push_back( idNow);
        }

        // Pick masses. Pick new channel if fail.
        if (!pickMasses()) continue;
        foundChannel = true;
        break;
      }

      // Failed to find acceptable decays.
      if (!foundChannel) {
        ostringstream osWarn;
        osWarn << "for id = " << id0;
        infoPtr->errorMsg("Error in ResonanceDecays::next:"
          " failed to find workable decay channel", osWarn.str());
        return false;
      }

      // Select colours in decay.
      if (!pickColours(iDec, process)) return false;

      // Select four-momenta in decay, boosted to lab frame.
      pProd.resize(0);
      pProd.push_back( decayer.p() );
      if (!pickKinematics()) return false;

      // Append decay products to the process event record. Set lifetimes.
      int iFirst = process.size();
        for (int i = 1; i <= mult; ++i) {
          process.append( idProd[i], 23, iDec, 0, 0, 0, cols[i], acols[i],
            pProd[i], mProd[i], m0);
        }
      int iLast = process.size() - 1;

      // Set decay vertex when this is displaced.
      if (process[iDec].hasVertex() || process[iDec].tau() > 0.) {
        Vec4 vDec = process[iDec].vDec();
        for (int i = iFirst; i <= iLast; ++i) process[i].vProd( vDec );
      }

      // Set lifetime of daughters.
      for (int i = iFirst; i <= iLast; ++i)
        process[i].tau( process[i].tau0() * rndmPtr->exp() );

      // Modify mother status and daughters.
      decayer.status(-22);
      decayer.daughters(iFirst, iLast);

    // End of loop over all entries.
    }
  } while (iDecNow == 0 && ++iDec < process.size());

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select masses of decay products.

bool ResonanceDecays::pickMasses() {

  // Arrays with properties of particles. Fill with dummy values for mother.
  vector<bool>   useBW;
  vector<double> m0BW, mMinBW, mMaxBW, widthBW;
  double mMother  = mProd[0];
  double m2Mother = mMother * mMother;
  useBW.push_back( false );
  m0BW.push_back( mMother );
  mMinBW.push_back( mMother );
  mMaxBW.push_back( mMother );
  widthBW.push_back( 0. );

  // Loop throught products for masses and widths. Set nominal mass.
  bool   useBWNow;
  double m0Now, mMinNow, mMaxNow, widthNow;
  for (int i = 1; i <= mult; ++i) {
    useBWNow  = particleDataPtr->useBreitWigner( idProd[i] );
    m0Now     = particleDataPtr->m0( idProd[i] );
    mMinNow   = particleDataPtr->m0Min( idProd[i] );
    mMaxNow   = particleDataPtr->m0Max( idProd[i] );
    if (useBWNow && mMaxNow < mMinNow) mMaxNow = mMother;
    widthNow  = particleDataPtr->mWidth( idProd[i] );
    useBW.push_back( useBWNow );
    m0BW.push_back( m0Now );
    mMinBW.push_back( mMinNow );
    mMaxBW.push_back( mMaxNow );
    widthBW.push_back( widthNow );
    mProd.push_back( m0Now );
  }

  // Find number of Breit-Wigners and summed (minimal) masses.
  int    nBW     = 0;
  double mSum    = 0.;
  double mSumMin = 0.;
  for (int i = 1; i <= mult; ++i) {
    if (useBW[i]) ++nBW;
    mSum        += max( m0BW[i], mMinBW[i]);
    mSumMin     += mMinBW[i];
  }

  // If sum of minimal masses above mother mass then give up.
  if (mSumMin + MSAFETY > mMother) return false;

  // If sum of masses below and no Breit-Wigners then done.
  if (mSum + 0.5 * MSAFETY < mMother && nBW == 0) return true;

  // Else if below then retry Breit-Wigners, with simple treshold.
  if (mSum + MSAFETY < mMother) {
    double wtMax = 2. * sqrtpos(1. - mSum*mSum / m2Mother);
    double wt;
    for (int iTryMasses = 0; iTryMasses <= NTRYMASSES; ++ iTryMasses) {
      if (iTryMasses == NTRYMASSES) return false;
      mSum = 0.;
      for (int i = 1; i <= mult; ++i) {
        if (useBW[i])  mProd[i] = particleDataPtr->mSel( idProd[i] );
        mSum += mProd[i];
      }
      wt = (mSum + 0.5 * MSAFETY < mMother)
         ? sqrtpos(1. - mSum*mSum / m2Mother) : 0.;
      if (wt > rndmPtr->flat() * wtMax) break;
    }
    return true;
  }

  // From now on some particles will have to be forced off shell.

  // Order Breit-Wigners in decreasing widths. Sum of other masses.
  vector<int> iBW;
  double mSum0 = 0.;
  for (int i = 1; i <= mult; ++i) {
    if (useBW[i]) iBW.push_back(i);
    else          mSum0 += mProd[i];
  }
  for (int i = 1; i < nBW; ++i) {
    for (int j = i - 1; j >= 0; --j) {
      if (widthBW[iBW[j+1]] > widthBW[iBW[j]]) swap (iBW[j+1], iBW[j]);
      else break;
    }
  }

  // Do all but broadest two in increasing-width order. Includes only one.
  if (nBW != 2) {
    int iMin = (nBW == 1) ? 0 : 2;
    for (int i = nBW - 1; i >= iMin; --i) {
      int iBWi = iBW[i];

      // Find allowed mass range of current resonance.
      double mMax    = mMother - mSum0 - MSAFETY;
      if (nBW  != 1) for (int j = 0; j < i; ++j) mMax -= mMinBW[iBW[j]];
      mMax           = min( mMaxBW[iBWi], mMax );
      double mMin    = min( mMinBW[iBWi], mMax - MSAFETY);
      if (mMin < 0.) return false;

      // Parameters for Breit-Wigner choice, with constrained mass range.
      double m2Nom   = pow2( m0BW[iBWi] );
      double m2Max   = mMax * mMax;
      double m2Min   = mMin * mMin;
      double mmWid   = m0BW[iBWi] * widthBW[iBWi];
      double atanMin = atan( (m2Min - m2Nom) / mmWid );
      double atanMax = atan( (m2Max - m2Nom) / mmWid );
      double atanDif = atanMax - atanMin;

      // Fail if too narrow mass range; e.g. out in tail of Breit-Wigner.
      if (atanDif < TINYBWRANGE) return false;

      // Retry mass according to Breit-Wigner, with simple threshold factor.
      double mr1     = mSum0*mSum0 / m2Mother;
      double mr2     = m2Min / m2Mother;
      double wtMax   = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );
      double m2Now, wt;
      for (int iTryMasses = 0; iTryMasses <= NTRYMASSES; ++ iTryMasses) {
        if (iTryMasses == NTRYMASSES) return false;
        m2Now = m2Nom + mmWid * tan(atanMin + rndmPtr->flat() * atanDif);
        mr2   = m2Now / m2Mother;
        wt    = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );
        if (wt > rndmPtr->flat() * wtMax) break;
      }

      // Prepare to iterate for more. Done for one Breit-Wigner.
      mProd[iBWi] = sqrt(m2Now);
      mSum0        += mProd[iBWi];
    }
    if (nBW == 1) return true;
  }

  // Left to do two broadest Breit-Wigners correlated, i.e. more realistic.
  int iBW1        = iBW[0];
  int iBW2        = iBW[1];
  int idMother    = abs(idProd[0]);
  int idDau1      = abs(idProd[iBW1]);
  int idDau2      = abs(idProd[iBW2]);

  // In some cases known phase-space behaviour; else simple beta factor.
  int psMode      = 1 ;
  if ( (idMother == 25 || idMother == 35) && idDau1 < 19
    && idDau2 == idDau1 ) psMode = 3;
  if ( (idMother == 25 || idMother == 35 )
    && (idDau1 == 23 || idDau1 == 24) && idDau2 == idDau1 ) psMode = 5;
  if ( idMother == 36
    && (idDau1 == 23 || idDau1 == 24) && idDau2 == idDau1 ) psMode = 6;

  // Find allowed mass ranges. Ensure that they are not closed.
  double mRem     = mMother - mSum0 - MSAFETY;
  double mMax1    = min( mMaxBW[iBW1], mRem - mMinBW[iBW2] );
  double mMin1    = min( mMinBW[iBW1], mMax1 - MSAFETY);
  double mMax2    = min( mMaxBW[iBW2], mRem - mMinBW[iBW1] );
  double mMin2    = min( mMinBW[iBW2], mMax2 - MSAFETY);

  // At least one range must extend below half remaining mass.
  if (mMin1 + mMin2 > mRem) return false;
  double mMid     = 0.5 * mRem;
  bool   hasMid1  = (mMin1 < mMid);
  bool   hasMid2  = (mMin2 < mMid);
  if (!hasMid1 && !hasMid2) return false;

  // Parameters for Breit-Wigner choice, with constrained mass range.
  double m2Nom1   = pow2( m0BW[iBW1] );
  double m2Max1   = mMax1 * mMax1;
  double m2Min1   = mMin1 * mMin1;
  double m2Mid1   = min( mMid * mMid, m2Max1);
  double mmWid1   = m0BW[iBW1] * widthBW[iBW1];
  double atanMin1 = atan( (m2Min1 - m2Nom1) / mmWid1 );
  double atanMax1 = atan( (m2Max1 - m2Nom1) / mmWid1 );
  double atanMid1 = (hasMid1) ? atan( (m2Mid1 - m2Nom1) / mmWid1 ) : 0.;
  double m2Nom2   = pow2( m0BW[iBW2] );
  double m2Max2   = mMax2 * mMax2;
  double m2Min2   = mMin2 * mMin2;
  double m2Mid2   = min( mMid * mMid, m2Max2);
  double mmWid2   = m0BW[iBW2] * widthBW[iBW2];
  double atanMin2 = atan( (m2Min2 - m2Nom2) / mmWid2 );
  double atanMax2 = atan( (m2Max2 - m2Nom2) / mmWid2 );
  double atanMid2 = (hasMid2) ? atan( (m2Mid2 - m2Nom2) / mmWid2 ) : 0.;

  // Relative weight to pick either below half remaining mass.
  double probLow1 = (hasMid1) ? 1. : 0.;
  if (hasMid1 && hasMid2) {
    double intLow1 = (atanMid1 - atanMin1) * (atanMax2 - atanMin2);
    double intLow2 = (atanMax1 - atanMin1) * (atanMid2 - atanMin2);
    probLow1 = intLow1 / (intLow1 + intLow2);
  }

  // Maximum matrix element times phase space weight.
  double m2Rem    = mRem * mRem;
  double mr1      = m2Min1 / m2Rem;
  double mr2      = m2Min2 / m2Rem;
  double psMax    = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );
  double wtMax   = 1.;
  if      (psMode == 1) wtMax = psMax;
  else if (psMode == 2) wtMax = psMax * psMax;
  else if (psMode == 3) wtMax = pow3(psMax);
  else if (psMode == 5) wtMax = psMax
    * (pow2(1. - mr1 - mr2) + 8. * mr1 * mr2);
  else if (psMode == 6) wtMax = pow3(psMax);

  // Retry mass according to Breit-Wigners, with simple threshold factor.
  double atanDif1, atanDif2, m2Now1, m2Now2, mNow1, mNow2, ps, wt;
  for (int iTryMasses = 0; iTryMasses <= NTRYMASSES; ++ iTryMasses) {
    if (iTryMasses == NTRYMASSES) return false;

    // Pick either below half remaining mass.
    bool pickLow1 = false;
    if (rndmPtr->flat() < probLow1) {
      atanDif1 = atanMid1 - atanMin1;
      atanDif2 = atanMax2 - atanMin2;
      pickLow1 = true;
    } else {
      atanDif1 = atanMax1 - atanMin1;
      atanDif2 = atanMid2 - atanMin2;
    }
    m2Now1 = m2Nom1 + mmWid1 * tan(atanMin1 + rndmPtr->flat() * atanDif1);
    m2Now2 = m2Nom2 + mmWid2 * tan(atanMin2 + rndmPtr->flat() * atanDif2);
    mNow1  = sqrt(m2Now1);
    mNow2  = sqrt(m2Now2);

    // Check that intended mass ordering is fulfilled.
    bool rejectRegion = (pickLow1) ? (mNow1 > mNow2) : (mNow2 > mNow1);
    if (rejectRegion) continue;

    // Threshold weight.
    mr1    = m2Now1 / m2Rem;
    mr2    = m2Now2 / m2Rem;
    wt     = 0.;
    if (mNow1 + mNow2 + MSAFETY < mMother) {
      ps   = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );
      wt   = 1.;
      if      (psMode == 1) wt = ps;
      else if (psMode == 2) wt = ps * ps;
      else if (psMode == 3) wt = pow3(ps);
      else if (psMode == 5) wt = ps
        * (pow2(1. - mr1 - mr2) + 8. * mr1 * mr2);
      else if (psMode == 6) wt = pow3(ps)*mr1*mr2;
    }
    if (wt > rndmPtr->flat() * wtMax) break;
  }
  mProd[iBW1] = mNow1;
  mProd[iBW2] = mNow2;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select colours of decay products.

bool ResonanceDecays::pickColours(int iDec, Event& process) {

  // Reset or create arrays with colour info.
  cols.resize(0);
  acols.resize(0);
  vector<int> iTriplet, iAtriplet, iOctet, iDipCol, iDipAcol;

  // Mother colours already known.
  int col0     = process[iDec].col();
  int acol0    = process[iDec].acol();
  int colType0 = process[iDec].colType();
  cols.push_back(  col0);
  acols.push_back(acol0);

  // Loop through all daughters.
  int colTypeNow;
  for (int i = 1; i <= mult; ++i) {
    // Daughter colours initially empty, so that all is set for singlet.
    cols.push_back(0);
    acols.push_back(0);
    // Find character (singlet, triplet, antitriplet, octet) of daughters.
    colTypeNow = particleDataPtr->colType( idProd[i] );
    if      (colTypeNow ==  0);
    else if (colTypeNow ==  1) iTriplet.push_back(i);
    else if (colTypeNow == -1) iAtriplet.push_back(i);
    else if (colTypeNow ==  2) iOctet.push_back(i);
    // Add two entries for sextets;
    else if (colTypeNow ==  3) {
      iTriplet.push_back(i);
      iTriplet.push_back(i);
    } else if (colTypeNow == -3) {
      iAtriplet.push_back(i);
      iAtriplet.push_back(i);
    } else {
      infoPtr->errorMsg("Error in ResonanceDecays::pickColours:"
        " unknown colour type encountered");
      return false;
    }
  }

  // Check excess of colours and anticolours in final over initial state.
  int nCol = iTriplet.size();
  if (colType0 == 1 || colType0 == 2) nCol -= 1;
  else if (colType0 == 3) nCol -= 2;
  int nAcol = iAtriplet.size();
  if (colType0 == -1 || colType0 == 2) nAcol -= 1;
  else if (colType0 == -3) nAcol -= 2;

  // If net creation of three colours then find junction kind:
  // mother is 1 = singlet, triplet, or sextet (no incoming RPV tags)
  //           3 = antitriplet, octet, or antisextet (acol0 = incoming RPV tag)
  //           5 = not applicable to decays (needs two incoming RPV tags)
  if (nCol - nAcol == 3) {
    int kindJun = (colType0 == 0 || colType0 == 1 || colType0 == 3) ? 1 : 3;

    // Set colours in three junction legs and store junction.
    int colJun[3];
    colJun[0] = (kindJun == 1) ? process.nextColTag() : acol0;
    colJun[1] = process.nextColTag();
    colJun[2] = process.nextColTag();
    process.appendJunction( kindJun, colJun[0], colJun[1], colJun[2]);

    // Loop over three legs. Remove an incoming anticolour on first leg.
    for (int leg = 0; leg < 3; ++leg) {
      if (leg == 0 && kindJun != 1) acol0 = 0;

      // Pick final-state triplets to carry these new colours.
      else {
        int pickT    = (iTriplet.size() == 1) ? 0
          : int( TINY + rndmPtr->flat() * (iTriplet.size() - TINY) );
        int iPickT   = iTriplet[pickT];
        cols[iPickT] = colJun[leg];

        // Remove matched triplet and store new colour dipole ends.
        iTriplet[pickT] = iTriplet.back();
        iTriplet.pop_back();
        iDipCol.push_back(iPickT);
        iDipAcol.push_back(0);
      }
    }

    // Update colour counter. Done with junction.
    nCol -= 3;
  }

  // If net creation of three anticolours then find antijunction kind:
  // mother is 2 = singlet, antitriplet, or antisextet (no incoming RPV tags)
  //           4 = triplet, octet, or sextet (col0 = incoming RPV tag)
  //           6 = not applicable to decays (needs two incoming RPV tags)
  if (nAcol - nCol == 3) {
    int kindJun = (colType0 == 0 || colType0 == -1 || colType0 == -3) ? 2 : 4;

    // Set anticolours in three antijunction legs and store antijunction.
    int acolJun[3];
    acolJun[0] = (kindJun == 2) ? process.nextColTag() : col0;
    acolJun[1] = process.nextColTag();
    acolJun[2] = process.nextColTag();
    process.appendJunction( kindJun, acolJun[0], acolJun[1], acolJun[2]);

    // Loop over three legs. Remove an incoming colour on first leg.
    for (int leg = 0; leg < 3; ++leg) {
      if (leg == 0 && kindJun != 2) col0 = 0;

      // Pick final-state antitriplets to carry these new anticolours.
      else {
        int pickA     = (iAtriplet.size() == 1) ? 0
          : int( TINY + rndmPtr->flat() * (iAtriplet.size() - TINY) );
        int iPickA    = iAtriplet[pickA];
        acols[iPickA] = acolJun[leg];

        // Remove matched antitriplet and store new colour dipole ends.
        iAtriplet[pickA] = iAtriplet.back();
        iAtriplet.pop_back();
        iDipCol.push_back(0);
        iDipAcol.push_back(iPickA);
      }
    }

    // Update anticolour counter. Done with antijunction.
    nAcol -= 3;
  }

  // If colours and anticolours do not match now then unphysical.
  if (nCol != nAcol) {
    infoPtr->errorMsg("Error in ResonanceDecays::pickColours:"
      " inconsistent colour tags");
    return false;
  }

  // Pick final-state triplet (if any) to carry initial colour.
  if (col0 > 0 && iTriplet.size() > 0) {
    int pickT    = (iTriplet.size() == 1) ? 0
      : int( TINY + rndmPtr->flat() * (iTriplet.size() - TINY) );
    int iPickT = iTriplet[pickT];
    cols[iPickT] = col0;

    // Remove matched triplet and store new colour dipole ends.
    col0 = 0;
    iTriplet[pickT] = iTriplet.back();
    iTriplet.pop_back();
    iDipCol.push_back(iPickT);
    iDipAcol.push_back(0);
  }

  // Pick final-state antitriplet (if any) to carry initial anticolour.
  if (acol0 > 0 && iAtriplet.size() > 0) {
    int pickA = (iAtriplet.size() == 1) ? 0
      : int( TINY + rndmPtr->flat() * (iAtriplet.size() - TINY) );
    int iPickA = iAtriplet[pickA];
    acols[iPickA] = acol0;

    // Remove matched antitriplet and store new colour dipole ends.
    acol0 = 0;
    iAtriplet[pickA] = iAtriplet.back();
    iAtriplet.pop_back();
    iDipCol.push_back(0);
    iDipAcol.push_back(iPickA);
  }

  // Sextets: second final-state triplet (if any)
  if (acol0 < 0 && iTriplet.size() > 0) {
    int pickT = (iTriplet.size() == 1) ? 0
      : int( TINY + rndmPtr->flat() * (iTriplet.size() - TINY) );
    int iPickT = iTriplet[pickT];
    cols[iPickT] = -acol0;

    // Remove matched antitriplet and store new colour dipole ends.
    acol0 = 0;
    iTriplet[pickT] = iTriplet.back();
    iTriplet.pop_back();
    iDipCol.push_back(iPickT);
    iDipAcol.push_back(0);
  }

  // Sextets: second final-state antitriplet (if any)
  if (col0 < 0 && iAtriplet.size() > 0) {
    int pickA    = (iAtriplet.size() == 1) ? 0
      : int( TINY + rndmPtr->flat() * (iAtriplet.size() - TINY) );
    int iPickA = iAtriplet[pickA];
    acols[iPickA] = -col0;

    // Remove matched triplet and store new colour dipole ends.
    col0 = 0;
    iAtriplet[pickA] = iAtriplet.back();
    iAtriplet.pop_back();
    iDipCol.push_back(0);
    iDipAcol.push_back(iPickA);
  }

  // Error checks that amount of leftover colours and anticolours match.
  if ( (iTriplet.size() != iAtriplet.size())
    || (col0 != 0 && acol0 == 0) || (col0 == 0 && acol0 != 0) ) {
    infoPtr->errorMsg("Error in ResonanceDecays::pickColours:"
      " inconsistent colour tags");
    return false;
  }

  // Match triplets to antitriplets in the final state.
  for (int pickT = 0; pickT < int(iTriplet.size()); ++pickT) {
    int iPickT = iTriplet[pickT];
    int pickA  = (iAtriplet.size() == 1) ? 0
      : int( TINY + rndmPtr->flat() * (iAtriplet.size() - TINY) );
    int iPickA = iAtriplet[pickA];

    // Connect pair with new colour tag.
    cols[iPickT]  = process.nextColTag();
    acols[iPickA] = cols[iPickT];

    // Remove matched antitriplet and store new colour dipole ends.
    iAtriplet[pickT] = iAtriplet.back();
    iAtriplet.pop_back();
    iDipCol.push_back(iPickT);
    iDipAcol.push_back(iPickA);
  }

  // If no octets are around then matching is done.
  if (col0 == 0 && acol0 == 0 && iOctet.size() == 0) return true;

  // If initial-state octet remains then store as (first!) new dipole.
  if (col0 != 0) {
    iDipCol.push_back(0);
    iDipAcol.push_back(0);
  }

  // Now attach all final-state octets at random to existing dipoles.
  for (int i = 0; i < int(iOctet.size()); ++i) {
    int iOct = iOctet[i];

    // If no dipole then start new one. (Happens for singlet -> octets.)
    if (iDipCol.size() == 0) {
      cols[iOct]  = process.nextColTag();
      acols[iOct] = cols[iOct] ;
      iDipCol.push_back(iOct);
      iDipAcol.push_back(iOct);
    }

    // Else attach to existing dipole picked at random.
    else {
      int pickDip = (iDipCol.size() == 1) ? 0
        : int( TINY + rndmPtr->flat() * (iDipCol.size() - TINY) );

      // Case with dipole in initial state: reattach existing colours.
      if (iDipCol[pickDip] == 0 && iDipAcol[pickDip] == 0) {
        cols[iOct]        = col0;
        acols[iOct]       = acol0;
        iDipAcol[pickDip] = iOct;
        iDipCol.push_back(iOct);
        iDipAcol.push_back(0);

      // Case with dipole from colour in initial state: also new colour.
      } else if (iDipAcol[pickDip] == 0) {
        int iPickCol      = iDipCol[pickDip];
        cols[iOct]        = cols[iPickCol];
        acols[iOct]       = process.nextColTag();
        cols[iPickCol]    = acols[iOct];
        iDipCol[pickDip]  = iOct;
        iDipCol.push_back(iPickCol);
        iDipAcol.push_back(iOct);

      // Remaining cases with dipole from anticolour in initial state
      // or dipole inside final state: also new colour.
      } else {
        int iPickAcol     = iDipAcol[pickDip];
        acols[iOct]       = acols[iPickAcol];
        cols[iOct]        = process.nextColTag();
        acols[iPickAcol]  = cols[iOct];
        iDipAcol[pickDip] = iOct;
        iDipCol.push_back(iOct);
        iDipAcol.push_back(iPickAcol);
      }
    }
  }

  // Must now have at least two dipoles (no 1 -> 8 or 8 -> 1).
  if (iDipCol.size() < 2) {
    infoPtr->errorMsg("Error in ResonanceDecays::pickColours:"
      " inconsistent colour tags");
    return false;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select decay products momenta isotropically in phase space.
// Process-dependent angular distributions may be imposed in SigmaProcess.

bool ResonanceDecays::pickKinematics() {

  // Description of two-body decays as simple special case.
  if (mult == 2) {

    // Masses.
    m0          = mProd[0];
    double m1   = mProd[1];
    double m2   = mProd[2];

    // Energies and absolute momentum in the rest frame.
    double e1   = 0.5 * (m0*m0 + m1*m1 - m2*m2) / m0;
    double e2   = 0.5 * (m0*m0 + m2*m2 - m1*m1) / m0;
    double pAbs = 0.5 * sqrtpos( (m0 - m1 - m2) * (m0 + m1 + m2)
      * (m0 + m1 - m2) * (m0 - m1 + m2) ) / m0;

    // Pick isotropic angles to give three-momentum.
    double cosTheta = 2. * rndmPtr->flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi      = 2. * M_PI * rndmPtr->flat();
    double pX       = pAbs * sinTheta * cos(phi);
    double pY       = pAbs * sinTheta * sin(phi);
    double pZ       = pAbs * cosTheta;

    // Fill four-momenta in mother rest frame and then boost to lab frame.
    pProd.push_back( Vec4(  pX,  pY,  pZ, e1) );
    pProd.push_back( Vec4( -pX, -pY, -pZ, e2) );
    pProd[1].bst( pProd[0] );
    pProd[2].bst( pProd[0] );

    // Done for two-body decay.
    return true;
  }

  // Description of three-body decays as semi-simple special case.
  if (mult == 3) {

    // Masses.
    m0             = mProd[0];
    double m1      = mProd[1];
    double m2      = mProd[2];
    double m3      = mProd[3];
    double mDiff   = m0 - (m1 + m2 + m3);

    // Kinematical limits for 2+3 mass. Maximum phase-space weight.
    double m23Min  = m2 + m3;
    double m23Max  = m0 - m1;
    double p1Max   = 0.5 * sqrtpos( (m0 - m1 - m23Min) * (m0 + m1 + m23Min)
      * (m0 + m1 - m23Min) * (m0 - m1 + m23Min) ) / m0;
    double p23Max  = 0.5 * sqrtpos( (m23Max - m2 - m3) * (m23Max + m2 + m3)
      * (m23Max + m2 - m3) * (m23Max - m2 + m3) ) / m23Max;
    double wtPSmax = 0.5 * p1Max * p23Max;

    // Pick an intermediate mass m23 flat in the allowed range.
    double wtPS, m23, p1Abs, p23Abs;
    do {
      m23 = m23Min + rndmPtr->flat() * mDiff;

      // Translate into relative momenta and find phase-space weight.
      p1Abs  = 0.5 * sqrtpos( (m0 - m1 - m23) * (m0 + m1 + m23)
        * (m0 + m1 - m23) * (m0 - m1 + m23) ) / m0;
      p23Abs = 0.5 * sqrtpos( (m23 - m2 - m3) * (m23 + m2 + m3)
        * (m23 + m2 - m3) * (m23 - m2 + m3) ) / m23;
      wtPS   = p1Abs * p23Abs;

    // If rejected, try again with new invariant masses.
    } while ( wtPS < rndmPtr->flat() * wtPSmax );

    // Set up m23 -> m2 + m3 isotropic in its rest frame.
    double cosTheta = 2. * rndmPtr->flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi      = 2. * M_PI * rndmPtr->flat();
    double pX       = p23Abs * sinTheta * cos(phi);
    double pY       = p23Abs * sinTheta * sin(phi);
    double pZ       = p23Abs * cosTheta;
    double e2       = sqrt( m2*m2 + p23Abs*p23Abs);
    double e3       = sqrt( m3*m3 + p23Abs*p23Abs);
    Vec4 p2(  pX,  pY,  pZ, e2);
    Vec4 p3( -pX, -pY, -pZ, e3);

    // Set up 0 -> 1 + 23 isotropic in its rest frame.
    cosTheta        = 2. * rndmPtr->flat() - 1.;
    sinTheta        = sqrt(1. - cosTheta*cosTheta);
    phi             = 2. * M_PI * rndmPtr->flat();
    pX              = p1Abs * sinTheta * cos(phi);
    pY              = p1Abs * sinTheta * sin(phi);
    pZ              = p1Abs * cosTheta;
    double e1       = sqrt( m1*m1 + p1Abs*p1Abs);
    double e23      = sqrt( m23*m23 + p1Abs*p1Abs);
    pProd.push_back( Vec4( pX, pY, pZ, e1) );

    // Boost 2 + 3 to the 0 rest frame and then boost to lab frame.
    Vec4 p23( -pX, -pY, -pZ, e23);
    p2.bst( p23 );
    p3.bst( p23 );
    pProd.push_back( p2 );
    pProd.push_back( p3 );
    pProd[1].bst( pProd[0] );
    pProd[2].bst( pProd[0] );
    pProd[3].bst( pProd[0] );

    // Done for three-body decay.
    return true;
  }

  // Do a multibody decay using the M-generator algorithm.

  // Mother and sum daughter masses.
  m0             = mProd[0];
  double mSum    = mProd[1];
  for (int i = 2; i <= mult; ++i) mSum += mProd[i];
  double mDiff   = m0 - mSum;

  // Begin setup of intermediate invariant masses.
  vector<double> mInv;
  for (int i = 0; i <= mult; ++i) mInv.push_back( mProd[i]);

  // Calculate the maximum weight in the decay.
  double wtPSmax = 1. / WTCORRECTION[mult];
  double mMax    = mDiff + mProd[mult];
  double mMin    = 0.;
  for (int i = mult - 1; i > 0; --i) {
    mMax        += mProd[i];
    mMin        += mProd[i+1];
    double mNow  = mProd[i];
    wtPSmax *= 0.5 * sqrtpos( (mMax - mMin - mNow) * (mMax + mMin + mNow)
    * (mMax + mMin - mNow) * (mMax - mMin + mNow) ) / mMax;
  }

  // Begin loop to find the set of intermediate invariant masses.
  vector<double> rndmOrd;
  double wtPS;
  do {
    wtPS  = 1.;

    // Find and order random numbers in descending order.
    rndmOrd.resize(0);
    rndmOrd.push_back(1.);
    for (int i = 1; i < mult - 1; ++i) {
      double rndm = rndmPtr->flat();
      rndmOrd.push_back(rndm);
      for (int j = i - 1; j > 0; --j) {
        if (rndm > rndmOrd[j]) swap( rndmOrd[j], rndmOrd[j+1] );
        else break;
      }
    }
    rndmOrd.push_back(0.);

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
  vector<Vec4> pInv;
  pInv.resize(mult + 1);
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
    pProd.push_back( Vec4( pX, pY, pZ, eHad) );
    pInv[i+1].p( -pX, -pY, -pZ, eInv);
  }
  pProd.push_back( pInv[mult] );

  // Boost decay products to the mother rest frame and on to lab frame.
  pInv[1] = pProd[0];
  for (int iFrame = mult - 1; iFrame > 0; --iFrame)
    for (int i = iFrame; i <= mult; ++i) pProd[i].bst(pInv[iFrame]);

  // Done for multibody decay.
  return true;

}

//==========================================================================

} // end namespace Pythia8
