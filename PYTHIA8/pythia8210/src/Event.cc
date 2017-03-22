// Event.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// Particle and Event classes, and some related global functions.

#include "Pythia8/Event.h"

namespace Pythia8 {

//==========================================================================

// Particle class.
// This class holds info on a particle in general.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Small number to avoid division by zero.
const double Particle::TINY = 1e-20;

//--------------------------------------------------------------------------

// Set pointer to the particle data species of the particle.

void Particle::setPDEPtr(ParticleDataEntry* pdePtrIn) {
  pdePtr = pdePtrIn; if (pdePtrIn != 0 || evtPtr == 0) return;
  pdePtr = (*evtPtr).particleDataPtr->particleDataEntryPtr( idSave);}

//--------------------------------------------------------------------------

// Functions for rapidity and pseudorapidity.

double Particle::y() const {
  double temp = log( ( pSave.e() + abs(pSave.pz()) ) / max( TINY, mT() ) );
  return (pSave.pz() > 0) ? temp : -temp;
}

double Particle::eta() const {
  double temp = log( ( pSave.pAbs() + abs(pSave.pz()) ) / max( TINY, pT() ) );
  return (pSave.pz() > 0) ? temp : -temp;
}

//--------------------------------------------------------------------------

// Method to find the index of the particle in the event record.

int Particle::index() const { if (evtPtr == 0) return -1;
  return (long(this) - long(&((*evtPtr)[0]))) / sizeof(Particle);
}

//--------------------------------------------------------------------------

// Trace the first and last copy of one and the same particle.

int Particle::iTopCopy() const {

  if (evtPtr == 0) return -1;
  int iUp = index();
  while ( iUp > 0 && (*evtPtr)[iUp].mother2() == (*evtPtr)[iUp].mother1()
    && (*evtPtr)[iUp].mother1() > 0) iUp = (*evtPtr)[iUp].mother1();
  return iUp;

}

int Particle::iBotCopy() const {

  if (evtPtr == 0) return -1;
  int iDn = index();
  while ( iDn > 0 && (*evtPtr)[iDn].daughter2() == (*evtPtr)[iDn].daughter1()
    && (*evtPtr)[iDn].daughter1() > 0) iDn = (*evtPtr)[iDn].daughter1();
  return iDn;

}

//--------------------------------------------------------------------------

// Trace the first and last copy of one and the same particle,
// also through shower branchings, making use of flavour matches.
// Stops tracing when this gives ambiguities.

int Particle::iTopCopyId( bool simplify) const {

  // Check that particle belongs to event record. Initialize.
  if (evtPtr == 0) return -1;
  int iUp = index();

  // Simple solution when only first and last mother are studied.
  if (simplify) for ( ; ; ) {
    int mother1up = (*evtPtr)[iUp].mother1();
    int id1up     = (mother1up > 0) ? (*evtPtr)[mother1up].id() : 0;
    int mother2up = (*evtPtr)[iUp].mother2();
    int id2up     = (mother2up > 0) ? (*evtPtr)[mother2up].id() : 0;
    if (mother2up != mother1up && id2up == id1up) return iUp;
    if (id1up != idSave && id2up != idSave) return iUp;
    iUp = (id1up == idSave) ?  mother1up : mother2up;
  }

  // Else full solution where all mothers are studied.
  for ( ; ; ) {
    int iUpTmp  = 0;
    vector<int> mothersTmp = (*evtPtr)[iUp].motherList();
    for (unsigned int i = 0; i < mothersTmp.size(); ++i)
    if ( (*evtPtr)[mothersTmp[i]].id() == idSave) {
      if (iUpTmp != 0) return iUp;
      iUpTmp = mothersTmp[i];
    }
    if (iUpTmp == 0) return iUp;
    iUp = iUpTmp;
  }

}

int Particle::iBotCopyId( bool simplify) const {

  // Check that particle belongs to event record. Initialize.
  if (evtPtr == 0) return -1;
  int iDn = index();

  // Simple solution when only first and last daughter are studied.
  if (simplify) for ( ; ; ) {
    int daughter1dn = (*evtPtr)[iDn].daughter1();
    int id1dn       = (daughter1dn > 0) ? (*evtPtr)[daughter1dn].id() : 0;
    int daughter2dn = (*evtPtr)[iDn].daughter2();
    int id2dn       = (daughter2dn > 0) ? (*evtPtr)[daughter2dn].id() : 0;
    if (daughter2dn != daughter1dn && id2dn == id1dn) return iDn;
    if (id1dn != idSave && id2dn != idSave) return iDn;
    iDn = (id1dn == idSave) ?  daughter1dn : daughter2dn;
  }

  // Else full solution where all daughters are studied.
  for ( ; ; ) {
    int iDnTmp  = 0;
    vector<int> daughtersTmp = (*evtPtr)[iDn].daughterList();
    for (unsigned int i = 0; i < daughtersTmp.size(); ++i)
    if ( (*evtPtr)[daughtersTmp[i]].id() == idSave) {
      if (iDnTmp != 0) return iDn;
      iDnTmp = daughtersTmp[i];
    }
    if (iDnTmp == 0) return iDn;
    iDn = iDnTmp;
  }

}

//--------------------------------------------------------------------------

// Find complete list of mothers.

vector<int> Particle::motherList() const {

  // Vector of all the mothers; created empty. Done if no event pointer.
  vector<int> motherVec;
  if (evtPtr == 0) return motherVec;

  // Special cases in the beginning, where the meaning of zero is unclear.
  int statusSaveAbs = abs(statusSave);
  if  (statusSaveAbs == 11 || statusSaveAbs == 12) ;
  else if (mother1Save == 0 && mother2Save == 0) motherVec.push_back(0);

  // One mother or a carbon copy.
  else if (mother2Save == 0 || mother2Save == mother1Save)
    motherVec.push_back(mother1Save);

  // A range of mothers from string fragmentation.
  else if ( (statusSaveAbs >  80 && statusSaveAbs <  90)
         || (statusSaveAbs > 100 && statusSaveAbs < 107) )
    for (int iRange = mother1Save; iRange <= mother2Save; ++iRange)
      motherVec.push_back(iRange);

  // Two separate mothers.
  else {
    motherVec.push_back( min(mother1Save, mother2Save) );
    motherVec.push_back( max(mother1Save, mother2Save) );
  }

  // Done.
  return motherVec;

}

//--------------------------------------------------------------------------

// Find complete list of daughters.

vector<int> Particle::daughterList() const {

  // Vector of all the daughters; created empty. Done if no event pointer.
  vector<int> daughterVec;
  if (evtPtr == 0) return daughterVec;

  // Simple cases: no or one daughter.
  if (daughter1Save == 0 && daughter2Save == 0) ;
  else if (daughter2Save == 0 || daughter2Save == daughter1Save)
    daughterVec.push_back(daughter1Save);

  // A range of daughters.
  else if (daughter2Save > daughter1Save)
    for (int iRange = daughter1Save; iRange <= daughter2Save; ++iRange)
      daughterVec.push_back(iRange);

  // Two separated daughters.
  else {
    daughterVec.push_back(daughter2Save);
    daughterVec.push_back(daughter1Save);
  }

  // Special case for two incoming beams: attach further
  // initiators and remnants that have beam as mother.
  if (abs(statusSave) == 12 || abs(statusSave) == 13) {
    int i = index();
    for (int iDau = i + 1; iDau < evtPtr->size(); ++iDau)
    if ((*evtPtr)[iDau].mother1() == i) {
      bool isIn = false;
      for (int iIn = 0; iIn < int(daughterVec.size()); ++iIn)
        if (iDau == daughterVec[iIn]) isIn = true;
      if (!isIn) daughterVec.push_back(iDau);
    }
  }

  // Done.
  return daughterVec;

}

//--------------------------------------------------------------------------

// Find complete list of sisters. Optionally traces up with iTopCopy
// and down with iBotCopy to give sisters at same level of evolution.

vector<int> Particle::sisterList(bool traceTopBot) const {

  // Vector of all the sisters; created empty.
  vector<int> sisterVec;
  if (evtPtr == 0 || abs(statusSave) == 11) return sisterVec;

  // Find all daughters of the mother.
  int iUp     = (traceTopBot) ? iTopCopy() : index();
  int iMother = (*evtPtr)[iUp].mother1();
  vector<int> daughterVec = (*evtPtr)[iMother].daughterList();

  // Copy all daughters, excepting the input particle itself.
  for (int j = 0; j < int(daughterVec.size()); ++j) {
    int iDau  = daughterVec[j];
    if (iDau != iUp) {
      int iDn = (traceTopBot) ? (*evtPtr)[iDau].iBotCopy() : iDau;
      sisterVec.push_back( iDn);
    }
  }

  // Done.
  return sisterVec;

}

//--------------------------------------------------------------------------

// Check whether a given particle is an arbitrarily-steps-removed
// mother to another. For the parton -> hadron transition, only
// first-rank hadrons are associated with the respective end quark.

bool Particle::isAncestor(int iAncestor) const {

  // Begin loop to trace upwards from the daughter.
  if (evtPtr == 0) return false;
  int iUp = index();
  int sizeNow = (*evtPtr).size();
  for ( ; ; ) {

    // If positive match then done.
    if (iUp == iAncestor) return true;

    // If out of range then failed to find match.
    if (iUp <= 0 || iUp > sizeNow) return false;

    // If unique mother then keep on moving up the chain.
    int mother1up = (*evtPtr)[iUp].mother1();
    int mother2up = (*evtPtr)[iUp].mother2();
    if (mother2up == mother1up || mother2up == 0) {iUp = mother1up; continue;}

    // If many mothers, except hadronization, then fail tracing.
    int statusUp = (*evtPtr)[iUp].statusAbs();
    if (statusUp < 81 || statusUp > 86) return false;

    // For hadronization step, fail if not first rank, else move up.
    if (statusUp == 82) {
      iUp = (iUp + 1 < sizeNow && (*evtPtr)[iUp + 1].mother1() == mother1up)
          ? mother1up : mother2up; continue;
    }
    if (statusUp == 83) {
      if ((*evtPtr)[iUp - 1].mother1() == mother1up) return false;
      iUp = mother1up; continue;
    }
    if (statusUp == 84) {
      if (iUp + 1 < sizeNow && (*evtPtr)[iUp + 1].mother1() == mother1up)
        return false;
      iUp = mother1up; continue;
    }

    // Fail for ministring -> one hadron and for junctions.
    return false;

  }
  // End of loop. Should never reach beyond here.
  return false;

}

//--------------------------------------------------------------------------

// Convert internal Pythia status codes to the HepMC status conventions.

int Particle::statusHepMC() const {

  // Positive codes are final particles. Status -12 are beam particles.
  if (statusSave > 0)    return 1;
  if (statusSave == -12) return 4;
  if (evtPtr == 0) return 0;

  // Hadrons, muons, taus that decay normally are status 2.
  if (isHadron() || abs(idSave) == 13 || abs(idSave) == 15) {
    // Particle should not decay into itself  (e.g. Bose-Einstein).
    if ( (*evtPtr)[daughter1Save].id() != idSave) {
      int statusDau = (*evtPtr)[daughter1Save].statusAbs();
      if (statusDau > 90 && statusDau < 95) return 2;
    }
  }

  // Other acceptable negative codes as their positive counterpart.
  if (statusSave <= -11 && statusSave >= -200) return -statusSave;

  // Unacceptable codes as 0.
  return 0;

}

//--------------------------------------------------------------------------

// Check if particle belonged to the final state at the Parton Level.

bool Particle::isFinalPartonLevel() const {

  if (index() >= evtPtr->savedPartonLevelSize) return false;
  if (statusSave > 0) return true;
  if (daughter1Save >= evtPtr->savedPartonLevelSize) return true;
  return false;

}

//--------------------------------------------------------------------------

// Recursively remove the decay products of a particle, update it to be
// undecayed, and update all mother/daughter indices to be correct.
// Warning: assumes that decay chains are nicely ordered.

bool Particle::undoDecay() {

  // Do not remove daughters of a parton, i.e. entry carrying colour.
  if (evtPtr == 0) return false;
  int i = index();
  if (i < 0 || i >= int((*evtPtr).size())) return false;
  if (colSave != 0 || acolSave != 0) return false;

  // Find range of daughters to remove.
  int dau1 = daughter1Save;
  if (dau1 == 0) return false;
  int dau2 = daughter2Save;
  if (dau2 == 0) dau2 = dau1;

  // Refuse if any of the daughters have other mothers.
  for (int j = dau1; j <= dau2; ++j) if ((*evtPtr)[j].mother1() != i
    || ((*evtPtr)[j].mother2() != i && (*evtPtr)[j].mother2() != 0) )
    return false;

  // Initialize range arrays for daughters and granddaughters.
  vector<int> dauBeg, dauEnd;
  dauBeg.push_back( dau1);
  dauEnd.push_back( dau2);

  // Begin recursive search through all decay chains.
  int iRange = 0;
  do {
    for (int j = dauBeg[iRange]; j <= dauEnd[iRange]; ++j)
    if ((*evtPtr)[j].status() < 0) {

      // Find new daughter range, if present.
      dau1 = (*evtPtr)[j].daughter1();
      if (dau1 == 0) return false;
      dau2 = (*evtPtr)[j].daughter2();
      if (dau2 == 0) dau2 = dau1;

      // Check if the range duplicates or contradicts existing ones.
      bool isNew = true;
      for (int iR = 0; iR < int(dauBeg.size()); ++iR) {
        if (dau1 == dauBeg[iR] && dau2 == dauEnd[iR]) isNew = false;
        else if (dau1 >= dauBeg[iR] && dau1 <= dauEnd[iR]) return false;
        else if (dau2 >= dauBeg[iR] && dau2 <= dauEnd[iR]) return false;
      }

      // Add new range where relevant. Keep ranges ordered.
      if (isNew) {
        dauBeg.push_back( dau1);
        dauEnd.push_back( dau2);
        for (int iR = int(dauBeg.size()) - 1; iR > 0; --iR) {
          if (dauBeg[iR] < dauBeg[iR - 1]) {
            swap( dauBeg[iR], dauBeg[iR - 1]);
            swap( dauEnd[iR], dauEnd[iR - 1]);
          } else break;
        }
      }

    // End of recursive search all decay chains.
    }
  } while (++iRange < int(dauBeg.size()));

  // Join adjacent ranges to reduce number of erase steps.
  if (int(dauBeg.size()) > 1) {
    int iRJ = 0;
    do {
      if (dauEnd[iRJ] + 1 == dauBeg[iRJ + 1]) {
        for (int iRB = iRJ + 1; iRB < int(dauBeg.size()) - 1; ++iRB)
          dauBeg[iRB] = dauBeg[iRB + 1];
        for (int iRE = iRJ; iRE < int(dauEnd.size()) - 1; ++iRE)
          dauEnd[iRE] = dauEnd[iRE + 1];
        dauBeg.pop_back();
        dauEnd.pop_back();
      } else ++iRJ;
    } while (iRJ < int(dauBeg.size()) - 1);
  }

  // Iterate over relevant ranges, from bottom up.
  for (int iR = int(dauBeg.size()) - 1; iR >= 0; --iR) {
    dau1 = dauBeg[iR];
    dau2 = dauEnd[iR];
    int nRem = dau2 - dau1 + 1;

    // Remove daughters in each range.
    (*evtPtr).remove( dau1, dau2);

    // Update subsequent history to account for removed indices.
    for (int j = 0; j < int((*evtPtr).size()); ++j) {
      if ((*evtPtr)[j].mother1() > dau2)
        (*evtPtr)[j].mother1( (*evtPtr)[j].mother1() - nRem );
      if ((*evtPtr)[j].mother2() > dau2)
        (*evtPtr)[j].mother2( (*evtPtr)[j].mother2() - nRem );
      if ((*evtPtr)[j].daughter1() > dau2)
        (*evtPtr)[j].daughter1( (*evtPtr)[j].daughter1() - nRem );
      if ((*evtPtr)[j].daughter2() > dau2)
        (*evtPtr)[j].daughter2( (*evtPtr)[j].daughter2() - nRem );
    }
  }

  // Update mother that has been undecayed.
  statusSave = abs(statusSave);
  daughter1Save = 0;
  daughter2Save = 0;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Particle name, with status but imposed maximum length -> may truncate.

string Particle::nameWithStatus(int maxLen) const {

  if (pdePtr == 0) return " ";
  string temp = (statusSave > 0) ? pdePtr->name(idSave)
    : "(" + pdePtr->name(idSave) + ")";
  while (int(temp.length()) > maxLen) {
    // Remove from end, excluding closing bracket and charge.
    int iRem = temp.find_last_not_of(")+-0");
    temp.erase(iRem, 1);
  }
  return temp;
}

//--------------------------------------------------------------------------

// Add offsets to mother and daughter pointers (must be non-negative).

void Particle::offsetHistory( int minMother, int addMother, int minDaughter,
  int addDaughter) {

  if (addMother < 0 || addDaughter < 0) return;
  if (  mother1Save > minMother  )   mother1Save += addMother;
  if (  mother2Save > minMother  )   mother2Save += addMother;
  if (daughter1Save > minDaughter) daughter1Save += addDaughter;
  if (daughter2Save > minDaughter) daughter2Save += addDaughter;

}

//--------------------------------------------------------------------------

// Add offsets to colour and anticolour (must be positive).

void Particle::offsetCol( int addCol) {

  if (addCol < 0) return;
  if ( colSave > 0)  colSave += addCol;
  if (acolSave > 0) acolSave += addCol;

}

//--------------------------------------------------------------------------

// Invariant mass of a pair and its square.
// (Not part of class proper, but tightly linked.)

double m(const Particle& pp1, const Particle& pp2) {
  double m2 = pow2(pp1.e() + pp2.e()) - pow2(pp1.px() + pp2.px())
     - pow2(pp1.py() + pp2.py()) - pow2(pp1.pz() + pp2.pz());
  return (m2 > 0. ? sqrt(m2) : 0.);
}

double m2(const Particle& pp1, const Particle& pp2) {
  double m2 = pow2(pp1.e() + pp2.e()) - pow2(pp1.px() + pp2.px())
     - pow2(pp1.py() + pp2.py()) - pow2(pp1.pz() + pp2.pz());
  return m2;
}

//==========================================================================

// Event class.
// This class holds info on the complete event record.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maxmimum number of mothers or daughter indices per line in listing.
const int Event::IPERLINE = 20;

//--------------------------------------------------------------------------

// Copy all information from one event record to another.

Event& Event::operator=( const Event& oldEvent) {

  // Do not copy if same.
  if (this != &oldEvent) {

    // Reset all current info in the event.
    clear();

    // Copy particle data table; needed for individual particles.
    particleDataPtr     = oldEvent.particleDataPtr;

    // Copy all the particles one by one.
    maxColTag = 100;
    for (int i = 0; i < oldEvent.size(); ++i) append( oldEvent[i] );

    // Copy all the junctions one by one.
    for (int i = 0; i < oldEvent.sizeJunction(); ++i)
      appendJunction( oldEvent.getJunction(i) );

    // Copy all other values.
    startColTag         = oldEvent.startColTag;
    maxColTag           = oldEvent.maxColTag;
    savedSize           = oldEvent.savedSize;
    savedJunctionSize   = oldEvent.savedJunctionSize;
    scaleSave           = oldEvent.scaleSave;
    scaleSecondSave     = oldEvent.scaleSecondSave;
    headerList          = oldEvent.headerList;

  // Done.
  }
  return *this;

}

//--------------------------------------------------------------------------

// Add a copy of an existing particle at the end of the event record;
// return index. Three cases, depending on sign of new status code:
// Positive: copy is viewed as daughter, status of original is negated.
// Negative: copy is viewed as mother, status of original is unchanged.
// Zero: the new is a perfect carbon copy (maybe to be changed later).

int Event::copy(int iCopy, int newStatus) {

  // Protect against attempt to copy negative entries (e.g., junction codes)
  // or entries beyond end of record.
  if (iCopy < 0 || iCopy >= size()) return -1;

  // Simple carbon copy.
  entry.push_back(entry[iCopy]);
  int iNew = entry.size() - 1;

  // Set up to make new daughter of old.
  if (newStatus > 0) {
    entry[iCopy].daughters(iNew,iNew);
    entry[iCopy].statusNeg();
    entry[iNew].mothers(iCopy, iCopy);
    entry[iNew].status(newStatus);

  // Set up to make new mother of old.
  } else if (newStatus < 0) {
    entry[iCopy].mothers(iNew,iNew);
    entry[iNew].daughters(iCopy, iCopy);
    entry[iNew].status(newStatus);
  }

  // Done.
  return iNew;

}

//--------------------------------------------------------------------------

// Print an event - special cases that rely on the general method.
// Not inline to make them directly callable in (some) debuggers.

void Event::list(int precision) const {
  list(false, false, cout, precision);
}

void Event::list(ostream& os, int precision) const {
  list(false, false, os, precision);
}

void Event::list(bool showScaleAndVertex, bool showMothersAndDaughters,
  int precision) const {
  list(showScaleAndVertex, showMothersAndDaughters, cout, precision);
}

//--------------------------------------------------------------------------

// Print an event.

void Event::list(bool showScaleAndVertex, bool showMothersAndDaughters,
  ostream& os, int precision) const {

  // Header.
  os << "\n --------  PYTHIA Event Listing  " << headerList << "----------"
     << "-------------------------------------------------\n \n    no    "
     << "    id   name            status     mothers   daughters     colou"
     << "rs      p_x        p_y        p_z         e          m \n";
  if (showScaleAndVertex)
    os << "                                    scale         pol          "
       << "                   xProd      yProd      zProd      tProd      "
       << " tau\n";

  // Precision. At high energy switch to scientific format for momenta.
  int prec = max( 3, precision);
  bool useFixed = (entry.empty() || entry[0].e() < 1e5);

  // Listing of complete event.
  Vec4 pSum;
  double chargeSum = 0.;
  for (int i = 0; i < int(entry.size()); ++i) {
    const Particle& pt = entry[i];

    // Basic line for a particle, always printed.
    os << setw(6) << i << setw(10) << pt.id() << "   " << left
       << setw(18) << pt.nameWithStatus(18) << right << setw(4)
       << pt.status() << setw(6) << pt.mother1() << setw(6)
       << pt.mother2() << setw(6) << pt.daughter1() << setw(6)
       << pt.daughter2() << setw(6) << pt.col() << setw(6) << pt.acol()
       << ( (useFixed) ? fixed : scientific ) << setprecision(prec)
       << setw(8+prec) << pt.px() << setw(8+prec) << pt.py()
       << setw(8+prec) << pt.pz() << setw(8+prec) << pt.e()
       << setw(8+prec) << pt.m() << "\n";

    // Optional extra line for scale value, polarization and production vertex.
    if (showScaleAndVertex)
      os << "                              " << setw(8+prec) << pt.scale()
         << " " << fixed << setprecision(prec) << setw(8+prec) << pt.pol()
         << "                        " << scientific << setprecision(prec)
         << setw(8+prec) << pt.xProd() << setw(8+prec) << pt.yProd()
         << setw(8+prec) << pt.zProd() << setw(8+prec) << pt.tProd()
         << setw(8+prec) << pt.tau() << "\n";

    // Optional extra line, giving a complete list of mothers and daughters.
    if (showMothersAndDaughters) {
      int linefill = 2;
      os << "                mothers:";
      vector<int> allMothers = pt.motherList();
      for (int j = 0; j < int(allMothers.size()); ++j) {
        os << " " <<  allMothers[j];
        if (++linefill == IPERLINE) {os << "\n                "; linefill = 0;}
      }
      os << ";   daughters:";
      vector<int> allDaughters = pt.daughterList();
      for (int j = 0; j < int(allDaughters.size()); ++j) {
        os << " " <<  allDaughters[j];
        if (++linefill == IPERLINE) {os << "\n                "; linefill = 0;}
      }
      if (linefill !=0) os << "\n";
    }

    // Extra blank separation line when each particle spans more than one line.
    if (showScaleAndVertex || showMothersAndDaughters) os << "\n";

    // Statistics on momentum and charge.
    if (entry[i].status() > 0) {
      pSum += entry[i].p();
      chargeSum += entry[i].charge();
    }
  }

  // Line with sum charge, momentum, energy and invariant mass.
  os << fixed << setprecision(3) << "                                   "
     << "Charge sum:" << setw(7) << chargeSum << "           Momentum sum:"
     << ( (useFixed) ? fixed : scientific ) << setprecision(prec)
     << setw(8+prec) << pSum.px() << setw(8+prec) << pSum.py()
     << setw(8+prec) << pSum.pz() << setw(8+prec) << pSum.e()
     << setw(8+prec) << pSum.mCalc() << "\n";

  // Listing finished.
  os << "\n --------  End PYTHIA Event Listing  ----------------------------"
     << "-------------------------------------------------------------------"
     << endl;
}

//--------------------------------------------------------------------------

// Erase junction stored in specified slot and move up the ones under.

void Event::eraseJunction(int i) {

  for (int j = i; j < int(junction.size()) - 1; ++j)
    junction[j] = junction[j + 1];
  junction.pop_back();

}

//--------------------------------------------------------------------------

// Print the junctions in an event.

void Event::listJunctions(ostream& os) const {

  // Header.
  os << "\n --------  PYTHIA Junction Listing  "
     << headerList.substr(0,30) << "\n \n    no  kind  col0  col1  col2 "
     << "endc0 endc1 endc2 stat0 stat1 stat2\n";

  // Loop through junctions in event and list them.
  for (int i = 0; i < sizeJunction(); ++i)
    os << setw(6) << i << setw(6) << kindJunction(i) << setw(6)
       << colJunction(i, 0) << setw(6) << colJunction(i, 1) << setw(6)
       << colJunction(i, 2) << setw(6) << endColJunction(i, 0) << setw(6)
       << endColJunction(i, 1) << setw(6) << endColJunction(i, 2) << setw(6)
       << statusJunction(i, 0) << setw(6) << statusJunction(i, 1) << setw(6)
       << statusJunction(i, 2) << "\n";

  // Alternative if no junctions. Listing finished.
  if (sizeJunction() == 0) os << "    no junctions present \n";
  os << "\n --------  End PYTHIA Junction Listing  --------------------"
     << "------" << endl;
}

//--------------------------------------------------------------------------

// Operator overloading allows to append one event to an existing one.

Event& Event::operator+=( const Event& addEvent) {

  // Find offsets. One less since won't copy line 0.
  int offsetIdx = entry.size() - 1;
  int offsetCol = maxColTag;

  // Add energy to zeroth line and calculate new invariant mass.
  entry[0].p( entry[0].p() + addEvent[0].p() );
  entry[0].m( entry[0].mCalc() );

  // Read out particles from line 1 (not 0) onwards.
  Particle temp;
  for (int i = 1; i < addEvent.size(); ++i) {
    temp = addEvent[i];

    // Add offset to nonzero mother, daughter and colour indices.
    if (temp.mother1() > 0) temp.mother1( temp.mother1() + offsetIdx );
    if (temp.mother2() > 0) temp.mother2( temp.mother2() + offsetIdx );
    if (temp.daughter1() > 0) temp.daughter1( temp.daughter1() + offsetIdx );
    if (temp.daughter2() > 0) temp.daughter2( temp.daughter2() + offsetIdx );
    if (temp.col() > 0) temp.col( temp.col() + offsetCol );
    if (temp.acol() > 0) temp.acol( temp.acol() + offsetCol );

    // Append particle to summed event.
    append( temp );
  }

  // Read out junctions one by one.
  Junction tempJ;
  int begCol, endCol;
  for (int i = 0; i < addEvent.sizeJunction(); ++i) {
    tempJ = addEvent.getJunction(i);

    // Add colour offsets to all three legs.
    for (int  j = 0; j < 3; ++j) {
      begCol = tempJ.col(j);
      endCol = tempJ.endCol(j);
      if (begCol > 0) begCol += offsetCol;
      if (endCol > 0) endCol += offsetCol;
      tempJ.cols( j, begCol, endCol);
    }

    // Append junction to summed event.
    appendJunction( tempJ );
  }

  // Set header that indicates character as sum of events.
  headerList = "(combination of several events)  -------";

  // Done.
  return *this;

}

//==========================================================================

} // end namespace Pythia8
