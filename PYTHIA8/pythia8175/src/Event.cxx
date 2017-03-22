// Event.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// Particle and Event classes, and some related global functions.

#include "Event.h"

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

void Event::list() const {
  list(false, false, cout);
}

void Event::list(ostream& os) const {
  list(false, false, os);
}

void Event::list(bool showScaleAndVertex, bool showMothersAndDaughters)
  const {
  list(showScaleAndVertex, showMothersAndDaughters, cout);
}

//--------------------------------------------------------------------------

// Print an event.

void Event::list(bool showScaleAndVertex, bool showMothersAndDaughters, 
  ostream& os) const {

  // Header.
  os << "\n --------  PYTHIA Event Listing  " << headerList << "----------"
     << "-------------------------------------------------\n \n    no    "
     << "    id   name            status     mothers   daughters     colou"
     << "rs      p_x        p_y        p_z         e          m \n"; 
  if (showScaleAndVertex) 
    os << "                                    scale         pol          "
       << "                   xProd      yProd      zProd      tProd      "
       << " tau\n";  

  // At high energy switch to scientific format for momenta.
  bool useFixed = (entry[0].e() < 1e5);

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
       << ( (useFixed) ? fixed : scientific ) << setprecision(3) 
       << setw(11) << pt.px() << setw(11) << pt.py() << setw(11) 
       << pt.pz() << setw(11) << pt.e() << setw(11) << pt.m() << "\n";

    // Optional extra line for scale value, polarization and production vertex.
    if (showScaleAndVertex) 
      os << "                              " << setw(11) << pt.scale() 
         << " " << fixed << setprecision(3) << setw(11) << pt.pol()   
         << "                        " << scientific << setprecision(3) 
         << setw(11) << pt.xProd() << setw(11) << pt.yProd() 
         << setw(11) << pt.zProd() << setw(11) << pt.tProd() 
         << setw(11) << pt.tau() << "\n";

    // Optional extra line, giving a complete list of mothers and daughters.
    if (showMothersAndDaughters) {
      int linefill = 2;
      os << "                mothers:";
      vector<int> allMothers = motherList(i);
      for (int j = 0; j < int(allMothers.size()); ++j) {
        os << " " <<  allMothers[j];
        if (++linefill == IPERLINE) {os << "\n                "; linefill = 0;}
      }
      os << ";   daughters:";
      vector<int> allDaughters = daughterList(i);
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
     << ( (useFixed) ? fixed : scientific ) << setprecision(3) 
     << setw(11) << pSum.px() << setw(11) << pSum.py() << setw(11) 
     << pSum.pz() << setw(11) << pSum.e() << setw(11) << pSum.mCalc() 
     << "\n";

  // Listing finished.
  os << "\n --------  End PYTHIA Event Listing  ----------------------------"
     << "-------------------------------------------------------------------"
     << endl;
}

//--------------------------------------------------------------------------

// Recursively remove the decay products of particle i, update it to be 
// undecayed, and update all mother/daughter indices to be correct.
// Warning: assumes that decay chains are nicely ordered. 

bool Event::undoDecay(int i) {

  // Do not remove daughters of a parton, i.e. entry carrying colour.
  if (i < 0 || i >= int(entry.size())) return false;
  if (entry[i].col() != 0 || entry[i].acol() != 0) return false;

  // Find range of daughters to remove.
  int dau1 = entry[i].daughter1();
  if (dau1 == 0) return false; 
  int dau2 = entry[i].daughter2();
  if (dau2 == 0) dau2 = dau1;

  // Refuse if any of the daughters have other mothers.
  for (int j = dau1; j <= dau2; ++j) if (entry[j].mother1() != i 
    || (entry[j].mother2() != i && entry[j].mother2() != 0) ) return false;

  // Initialize range arrays for daughters and granddaughters.
  vector<int> dauBeg, dauEnd;
  dauBeg.push_back( dau1);
  dauEnd.push_back( dau2); 

  // Begin recursive search through all decay chains.
  int iRange = 0;
  do {
    for (int j = dauBeg[iRange]; j <= dauEnd[iRange]; ++j) 
    if (entry[j].status() < 0) {
      
      // Find new daughter range, if present.
      dau1 = entry[j].daughter1();
      if (dau1 == 0) return false; 
      dau2 = entry[j].daughter2();
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

  // Iterate over relevant ranges, from bottom up.
  for (int iR = int(dauBeg.size()) - 1; iR >= 0; --iR) {
    dau1 = dauBeg[iR];
    dau2 = dauEnd[iR];
    int nRem = dau2 - dau1 + 1;

    // Remove daughters in each range.
    entry.erase( entry.begin() + dau1, entry.begin() + dau2 + 1);

    // Update subsequent history to account for removed indices.
    for (int j = 0; j < int(entry.size()); ++j) {
      if (entry[j].mother1() > dau2)
        entry[j].mother1( entry[j].mother1() - nRem );  
      if (entry[j].mother2() > dau2)
        entry[j].mother2( entry[j].mother2() - nRem );  
      if (entry[j].daughter1() > dau2) 
        entry[j].daughter1( entry[j].daughter1() - nRem );  
      if (entry[j].daughter2() > dau2) 
        entry[j].daughter2( entry[j].daughter2() - nRem ); 
    }
  } 

  // Update mother that has been undecayed.
  entry[i].statusPos();
  entry[i].daughters();

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Find complete list of mothers.

vector<int> Event::motherList(int i) const {

  // Vector of all the mothers; created empty.
  vector<int> mothers;

  // Read out the two official mother indices and status code.
  int mother1   = entry[i].mother1();
  int mother2   = entry[i].mother2();
  int statusAbs = entry[i].statusAbs();

  // Special cases in the beginning, where the meaning of zero is unclear.
  if  (statusAbs == 11 || statusAbs == 12) ;
  else if (mother1 == 0 && mother2 == 0) mothers.push_back(0);
    
  // One mother or a carbon copy 
  else if (mother2 == 0 || mother2 == mother1) mothers.push_back(mother1); 

  // A range of mothers from string fragmentation.
  else if ( (statusAbs >  80 && statusAbs <  90)  
         || (statusAbs > 100 && statusAbs < 107) )  
    for (int iRange = mother1; iRange <= mother2; ++iRange) 
      mothers.push_back(iRange); 

  // Two separate mothers.
  else {
    mothers.push_back( min(mother1, mother2) ); 
    mothers.push_back( max(mother1, mother2) );
  }

  // Done.       
  return mothers;

}

//--------------------------------------------------------------------------

// Find complete list of daughters.

vector<int> Event::daughterList(int i) const {

  // Vector of all the daughters; created empty.
  vector<int> daughters;

  // Read out the two official daughter indices.
  int daughter1 = entry[i].daughter1();
  int daughter2 = entry[i].daughter2();

  // Simple cases: no or one daughter.
  if (daughter1 == 0 && daughter2 == 0) ;
  else if (daughter2 == 0 || daughter2 == daughter1) 
    daughters.push_back(daughter1);

  // A range of daughters.
  else if (daughter2 > daughter1)
    for (int iRange = daughter1; iRange <= daughter2; ++iRange) 
      daughters.push_back(iRange); 

  // Two separated daughters.
  else {
    daughters.push_back(daughter2); 
    daughters.push_back(daughter1);
  }

  // Special case for two incoming beams: attach further 
  // initiators and remnants that have beam as mother.
  if (entry[i].statusAbs() == 12 || entry[i].statusAbs() == 13)
  for (int iDau = i + 1; iDau < size(); ++iDau) 
  if (entry[iDau].mother1() == i) {
    bool isIn = false;
    for (int iIn = 0; iIn < int(daughters.size()); ++iIn) 
      if (iDau == daughters[iIn]) isIn = true;   
    if (!isIn) daughters.push_back(iDau); 
  }
    
  // Done.
  return daughters;

}

//--------------------------------------------------------------------------

// Convert internal Pythia status codes to the HepMC status conventions.

int Event::statusHepMC(int i) const {
  
  // Positive codes are final particles. Status -12 are beam particles.
  int statusNow     = entry[i].status();
  if (statusNow > 0)    return 1;
  if (statusNow == -12) return 4;  

  // Hadrons, muons, taus that decay normally are status 2.
  int idNow         = entry[i].id();
  if (entry[i].isHadron() || abs(idNow) == 13 || abs(idNow) == 15) {
    int iDau        = entry[i].daughter1();
    // Particle should not decay into itself  (e.g. Bose-Einstein).
    if ( entry[iDau].id() != idNow) {
      int statusDau = entry[ iDau ].statusAbs();
      if (statusDau > 90 && statusDau < 95) return 2;
    }
  }

  // Other acceptable negative codes as their positive counterpart.
  if (statusNow <= -11 && statusNow >= -200) return -statusNow;

  // Unacceptable codes as 0.
  return 0;

}

//--------------------------------------------------------------------------

// Trace the first and last copy of one and the same particle.

int Event::iTopCopy( int i) const {

  int iUp = i;
  while ( iUp > 0 && entry[iUp].mother2() == entry[iUp].mother1()
    && entry[iUp].mother1() > 0) iUp = entry[iUp].mother1();  
  return iUp;

}

int Event::iBotCopy( int i) const {

  int iDn = i;
  while ( iDn > 0 && entry[iDn].daughter2() == entry[iDn].daughter1()
    && entry[iDn].daughter1() > 0) iDn = entry[iDn].daughter1();  
  return iDn;

}

//--------------------------------------------------------------------------

// Trace the first and last copy of one and the same particle,
// also through shower branchings, making use of flavour matches.
// Stops tracing when this gives ambiguities.

int Event::iTopCopyId( int i) const {

  int id = entry[i].id();
  int iUp = i;
  for ( ; ; ) {
    int mother1 = entry[iUp].mother1();
    int id1     = (mother1 > 0) ? entry[mother1].id() : 0;
    int mother2 = entry[iUp].mother2();
    int id2     = (mother2 > 0) ? entry[mother2].id() : 0;
    if (mother2 != mother1 && id2 == id1) break;
    if (id1 == id) { 
      iUp = mother1; 
      continue;
    }
    if (id2 == id) {
      iUp = mother2; 
      continue;
    }
    break;
  } 
  return iUp;

}

int Event::iBotCopyId( int i) const {

  int id = entry[i].id();
  int iDn = i;
  for ( ; ; ) {
    int daughter1 = entry[iDn].daughter1();
    int id1       = (daughter1 > 0) ? entry[daughter1].id() : 0;
    int daughter2 = entry[iDn].daughter2();
    int id2       = (daughter2 > 0) ? entry[daughter2].id() : 0;
    if (daughter2 != daughter1 && id2 == id1) break;
    if (id1 == id) { 
      iDn = daughter1; 
      continue;
    }
    if (id2 == id) {
      iDn = daughter2; 
      continue;
    }
    break;
  } 
  return iDn;

}

//--------------------------------------------------------------------------

// Find complete list of sisters.

vector<int> Event::sisterList(int i) const {

  // Vector of all the sisters; created empty.
  vector<int> sisters;
  if (entry[i].statusAbs() == 11) return sisters;

  // Find mother and all its daughters.
  int iMother = entry[i].mother1();
  vector<int> daughters = daughterList(iMother);

  // Copy all daughters, excepting the input particle itself.
  for (int j = 0; j < int(daughters.size()); ++j) 
  if (daughters[j] != i) sisters.push_back( daughters[j] );

  // Done.       
  return sisters;

}

//--------------------------------------------------------------------------

// Find complete list of sisters. Traces up with iTopCopy and
// down with iBotCopy to give sisters at same level of evolution.
// Should this not give any final particles, the search is widened.

vector<int> Event::sisterListTopBot(int i, bool widenSearch) const {

  // Vector of all the sisters; created empty.
  vector<int> sisters;
  if (entry[i].statusAbs() == 11) return sisters;

  // Trace up to first copy of current particle.
  int iUp = iTopCopy(i);

  // Find mother and all its daughters.
  int iMother = entry[iUp].mother1();
  vector<int> daughters = daughterList(iMother);

  // Trace all daughters down, excepting the input particle itself.
  for (int jD = 0; jD < int(daughters.size()); ++jD) 
  if (daughters[jD] != iUp) 
    sisters.push_back( iBotCopy( daughters[jD] ) );

  // Prune any non-final particles from list.
  int jP = 0;
  while (jP < int(sisters.size())) {
    if (entry[sisters[jP]].status() > 0) ++jP;
    else {
      sisters[jP] = sisters.back();
      sisters.pop_back();
    }
  } 

  // If empty list then restore immediate daughters.
  if (sisters.size() == 0 && widenSearch) { 
    for (int jR = 0; jR < int(daughters.size()); ++jR) 
    if (daughters[jR] != iUp) 
      sisters.push_back( iBotCopy( daughters[jR] ) );
    
    // Then trace all daughters, not only bottom copy.
    for (int jT = 0; jT < int(sisters.size()); ++jT) {
      daughters = daughterList( sisters[jT] );
      for (int k = 0; k < int(daughters.size()); ++k)    
        sisters.push_back( daughters[k] );
    }
  
    // And then prune any non-final particles from list.
    int jN = 0;
    while (jN < int(sisters.size())) {
      if (entry[sisters[jN]].status() > 0) ++jN;
      else {
        sisters[jN] = sisters.back();
        sisters.pop_back();
      }
    } 
  }

  // Done.       
  return sisters;

}

//--------------------------------------------------------------------------

// Check whether a given particle is an arbitrarily-steps-removed
// mother to another. For the parton -> hadron transition, only 
// first-rank hadrons are associated with the respective end quark.

bool Event::isAncestor(int i, int iAncestor) const {

  // Begin loop to trace upwards from the daughter.
  int iUp = i;
  for ( ; ; ) {

    // If positive match then done.
    if (iUp == iAncestor) return true;

    // If out of range then failed to find match.
    if (iUp <= 0 || iUp > size()) return false;

    // If unique mother then keep on moving up the chain.
    int mother1 = entry[iUp].mother1();
    int mother2 = entry[iUp].mother2();
    if (mother2 == mother1 || mother2 == 0) {iUp = mother1; continue;}

    // If many mothers, except hadronization, then fail tracing.
    int status = entry[iUp].statusAbs();
    if (status < 81 || status > 86) return false;

    // For hadronization step, fail if not first rank, else move up.
    if (status == 82) {
      iUp = (iUp + 1 < size() && entry[iUp + 1].mother1() == mother1) 
          ? mother1 : mother2; continue;
    } 
    if (status == 83) {
      if (entry[iUp - 1].mother1() == mother1) return false;
      iUp = mother1; continue;
    }
    if (status == 84) {
      if (iUp + 1 < size() && entry[iUp + 1].mother1() == mother1) 
        return false;
      iUp = mother1; continue;
    }

    // Fail for ministring -> one hadron and for junctions. 
    return false;

  }
  // End of loop. Should never reach beyond here.
  return false;

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
