// ColosurReconnection.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// ColourReconnection class.

#include "Pythia8/ColourReconnection.h"

namespace Pythia8 {

//==========================================================================

// The BeamDipole class is purely internal to reconnectMPIs.

class BeamDipole {

public:

  // Constructor.
  BeamDipole( int colIn = 0, int iColIn = 0, int iAcolIn = 0)
    : col(colIn), iCol(iColIn), iAcol(iAcolIn) {}

  // Members.
  int    col, iCol, iAcol;
  double p1p2;

};

//==========================================================================

// The ColourDipole class.

//--------------------------------------------------------------------------

// Printing function, inteded for debugging.

void ColourDipole::print() {

  cout << setw(10) << this << setw(6) << col << setw(3) << colReconnection
       << setw(6) << iCol << setw(5) << iAcol << setw(6) << iColLeg << setw(5)
       << iAcolLeg << setw(6) << isJun << setw(5) << isAntiJun  << setw(10)
       << p1p2 << " colDips: ";
  for (int i = 0;i < int(colDips.size());++i)
    cout << setw(10) << colDips[i];
  cout  <<  " acolDips: ";
  for (int i = 0;i < int(acolDips.size());++i)
    cout << setw(10) << acolDips[i];
  cout << setw(3) << isActive << endl;

}

//==========================================================================

// The InfoGluonMove class is purely internal to reconnectMove.

class InfoGluonMove{

public:

  // Constructors.
  InfoGluonMove(int i1in, int col1in, int acol1in, int iCol1in, int iAcol1in,
    int col2in, int iCol2in, int iAcol2in, double lambdaRefIn,
    double dLambdaIn) : i1(i1in), i2(0), col1(col1in), acol1(acol1in),
    iCol1(iCol1in), iAcol1(iAcol1in), col2(col2in), iCol2(iCol2in),
    iAcol2(iAcol2in), lambdaRef(lambdaRefIn), dLambda(dLambdaIn) {}
  InfoGluonMove(int i1in, int i2in, int iCol1in, int iAcol1in, int iCol2in,
    int iAcol2in, int dLambdaIn) : i1(i1in), i2(i2in), col1(0), acol1(0),
    iCol1(iCol1in), iAcol1(iAcol1in), col2(0), iCol2(iCol2in),
    iAcol2(iAcol2in), lambdaRef(0.), dLambda(dLambdaIn) {}

  // Members.
  int i1, i2, col1, acol1, iCol1, iAcol1, col2, iCol2, iAcol2;
  double lambdaRef, dLambda;

};

//==========================================================================

// The ColourJunction class.

//--------------------------------------------------------------------------

// Printing function, inteded for debugging.

void ColourJunction::print() {

  cout << setw(6) << kind() << setw(6)
       << col(0) << setw(6) << col(1) << setw(6) << col(2) << setw(6)
       << endCol(0) << setw(6) << endCol(1) << setw(6) << endCol(2) << setw(6)
       << status(0) << setw(6) << status(1) << setw(6) << status(2) << setw(10)
       << dips[0] << setw(10) << dips[1] << setw(10) << dips[2] << setw(10)
       << "\n";
  cout << "     " << setw(10) << dipsOrig[0] << setw(10) << dipsOrig[1]
       << setw(10) << dipsOrig[2] << endl;

}

//==========================================================================

// The ColourParticle class.

//--------------------------------------------------------------------------

// Printing function, inteded for debugging.

void ColourParticle::list() {

  const Particle& pt = (*this);

  // Basic line for a particle, always printed.
  cout << setw(10) << pt.id() << "   " << left
       << setw(18) << pt.nameWithStatus(18) << right << setw(4)
       << pt.status() << setw(6) << pt.mother1() << setw(6)
       << pt.mother2() << setw(6) << pt.daughter1() << setw(6)
       << pt.daughter2() << setw(6) << pt.col() << setw(6) << pt.acol()
       << setprecision(3)
       << setw(11) << pt.px() << setw(11) << pt.py() << setw(11)
       << pt.pz() << setw(11) << pt.e() << setw(11) << pt.m() << "\n";

}

//--------------------------------------------------------------------------

// Printing function, inteded for debugging.

void ColourParticle::listActiveDips() {

  cout << "active dips: " << endl;
  for (int i = 0; i < int(activeDips.size()); ++i)
    activeDips[i]->print();

}

//--------------------------------------------------------------------------

// Printing function, inteded for debugging.

void ColourParticle::print() {

  cout << "---   Particle   ---" << endl;
  for (int i = 0; i < int(dips.size()); ++i) {
    cout << "(" <<colEndIncluded[i] << ") ";
    for (int j = 0; j < int(dips[i].size()); ++j) {
      cout << dips[i][j]->iCol << " (" << dips[i][j]->col << ") ";
      if (j == int(dips[i].size() - 1))
        cout << dips[i][j]->iAcol << " (" << acolEndIncluded[i] << ")" << endl;
    }
  }

}

//==========================================================================

// The ColourReconnection class.

//--------------------------------------------------------------------------

// Simple comparison function for sort.

bool cmpTrials(TrialReconnection j1, TrialReconnection j2) {
    return (j1.lambdaDiff < j2.lambdaDiff);
}

//--------------------------------------------------------------------------

// Initialization.

bool ColourReconnection::init( Info* infoPtrIn, Settings& settings,
  Rndm* rndmPtrIn,  BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  PartonSystems* partonSystemsPtrIn) {

  // Save pointers.
  infoPtr             = infoPtrIn;
  rndmPtr             = rndmPtrIn;
  beamAPtr            = beamAPtrIn;
  beamBPtr            = beamBPtrIn;
  partonSystemsPtr    = partonSystemsPtrIn;

  // Total and squared CM energy at nominal energy.
  eCM                 = infoPtr->eCM();
  sCM                 = eCM * eCM;

  // Choice of reconnection model.
  reconnectMode       = settings.mode("ColourReconnection:mode");

  // pT0 scale of MPI; used in the MPI-based reconnection model.
  pT0Ref              = settings.parm("MultipartonInteractions:pT0Ref");
  ecmRef              = settings.parm("MultipartonInteractions:ecmRef");
  ecmPow              = settings.parm("MultipartonInteractions:ecmPow");
  pT0                 = pT0Ref * pow(eCM / ecmRef, ecmPow);

  // Parameter of the MPI-based reconnection model.
  reconnectRange      = settings.parm("ColourReconnection:range");
  pT20Rec             = pow2(reconnectRange * pT0);

  // Parameters of the new reconnection model.
  m0                  = settings.parm("ColourReconnection:m0");
  m0sqr               = pow2(m0);
  allowJunctions      = settings.flag("ColourReconnection:allowJunctions");
  nReconCols          = settings.mode("ColourReconnection:nColours");
  sameNeighbourCol  = settings.flag("ColourReconnection:sameNeighbourColours");
  minimumGain         = 1e-10;
  minimumGainJun      = settings.parm("ColourReconnection:minimumGainJun");

  // Parameters of gluon-move model.
  m2Lambda            = settings.parm("ColourReconnection:m2Lambda");
  fracGluon           = settings.parm("ColourReconnection:fracGluon");
  dLambdaCut          = settings.parm("ColourReconnection:dLambdaCut");
  flipMode            = settings.mode("ColourReconnection:flipMode");

  // Initialize StringLength class.
  stringLength.init(infoPtr, settings);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do colour reconnection for current event.

bool ColourReconnection::next( Event& event, int iFirst) {

  // MPI-based reconnection model.
  if (reconnectMode == 0) return reconnectMPIs(event, iFirst);

  // New reconnection model.
  else if (reconnectMode == 1) return nextNew(event, iFirst);

  // Gluon-move model.
  else if (reconnectMode == 2) return reconnectMove(event, iFirst);

  // Undefined.
  else {
    infoPtr->errorMsg("Warning in ColourReconnection::next: "
                      "Colour reconnecion mode not found");
    return true;
  }

}

//--------------------------------------------------------------------------

// Do new colour reconnection for current event.

bool ColourReconnection::nextNew( Event& event, int iFirst) {

  // Clear old records.
  while (!dipoles.empty()) {
    delete dipoles.back();
    dipoles.pop_back();
  }
  particles.clear();
  junctions.clear();
  junTrials.clear();
  dipTrials.clear();

  // Setup dipoles and make pseudo particles.
  setupDipoles(event, iFirst);
  makeAllPseudoParticles(event, iFirst);

  // Setup all dipole reconnections.
  // Split dipoles into the 9 different "colours".
  vector<vector<int> > iDips;
  iDips.resize(nReconCols);
  for (int i = 0; i < int(iDips.size()); ++i)
    iDips[i] = vector<int>();

  for (int i = 0; i < int(dipoles.size()); ++i)
    if (dipoles[i]->isActive)
      iDips[dipoles[i]->colReconnection].push_back(i);

  // Loop over each colour individually.
  for (int i = 0;i < int(iDips.size()); ++i)
    for (int j = 0; j < int(iDips[i].size()); ++j)
      for (int k = j + 1; k < int(iDips[i].size()); ++k)
        singleReconnection(dipoles[iDips[i][j]], dipoles[iDips[i][k]]);

  // Start outer loop over reconnections.
  for (int iOuterLoop = 0; iOuterLoop < 20; ++iOuterLoop) {
    bool finished = true;

    // Do inner loop for string reconnections.
    for (int iTry = 0;iTry < 500 && dipTrials.size(); ++iTry) {

      // Store all dipoles connected to the chosen dipole.
      usedDipoles.clear();
      storeUsedDips(dipTrials.back());

      // Do the reconnection.
      doDipoleTrial(dipTrials.back());

      // Sort the used dipoles and remove copies of the same.
      sort(usedDipoles.begin(), usedDipoles.end());
      for (int i = 0;i < int(usedDipoles.size() - 1); ++i)
        if (usedDipoles[i] == usedDipoles[i + 1]) {
          usedDipoles.erase(usedDipoles.begin() + i);
          i--;
        }

      // Updating the dipole trials.
      updateDipoleTrials();

    }

    // Loop over list of dipoles to try and form junction structures.
    if (allowJunctions) {

      // Split dipoles into three categories.
      iDips.clear();
      iDips.resize(3);
      for (int i = 0; i < int(iDips.size()); ++i)
        iDips[i] = vector<int>();

      for (int i = 0; i < int(dipoles.size()); ++i)
        if (dipoles[i]->isActive)
          iDips[dipoles[i]->colReconnection % 3].push_back(i);

      // Loop over different "colours" (now only three different groups).
      for (int i = 0;i < int(iDips.size()); ++i)
        for (int j = 0; j < int(iDips[i].size()); ++j)
          for (int k = j + 1; k < int(iDips[i].size()); ++k)
            singleJunction(dipoles[iDips[i][j]], dipoles[iDips[i][k]]);


      // Loop over different "colours" (now only three different groups).
      for (int i = 0;i < int(iDips.size()); ++i)
        for (int j = 0; j < int(iDips[i].size()); ++j)
          for (int k = j + 1; k < int(iDips[i].size()); ++k)
            for (int l = k + 1; l < int(iDips[i].size()); ++l)
              singleJunction(dipoles[iDips[i][j]], dipoles[iDips[i][k]],
                                dipoles[iDips[i][l]]);                  

      for (int iTry = 0; iTry < 100 && junTrials.size(); ++iTry) {

        // Find all dipoles connected to the reconnection.
        usedDipoles.clear();
        storeUsedDips(junTrials.back());

        // Do the reconnection.
        doJunctionTrial(event, junTrials.back());

        // Sort the used dipoles and remove copies of the same.
        sort(usedDipoles.begin(), usedDipoles.end());
        for (int i = 0;i < int(usedDipoles.size() - 1); ++i)
          if (usedDipoles[i] == usedDipoles[i + 1]) {
            usedDipoles.erase(usedDipoles.begin() + i);
            i--;
          }

        // Update lists.
        updateJunctionTrials();
        updateDipoleTrials();

        finished = false;
      }
    }

    // If no junctions were made, the overall loop is finished.
    if (finished)
      break;
  }

  updateEvent(event, iFirst);

  // Done.
  return true;
}


//--------------------------------------------------------------------------

// Setup initial guess on dipoles, here all colours are assumed
// to be found in the final state.

void ColourReconnection::setupDipoles( Event& event, int iFirst) {

  // Make vectors needed for storage of chains.
  vector< vector<int > > chains;
  vector<bool> isJun;
  vector<bool> isAntiJun;
  vector<bool> isGluonLoop;
  vector<bool> inChain(event.size(),false);

  // Find all quarks and follow untill no more colour.
  for (int i = iFirst; i < event.size(); ++i) {
    if ( (event[i].isFinal() && !inChain[i]
      &&  event[i].isQuark() && event[i].id() > 0)
      || (event[i].isFinal() && !inChain[i]
      &&  event[i].isDiquark() && event[i].id() < 0) ) {
      int curCol = event[i].col();
      inChain[i] = true;
      vector<int> chain;
      chain.push_back(i);
      isAntiJun.push_back(false);
      isJun.push_back(false);
      isGluonLoop.push_back(false);
      for (int iSteps = 0; curCol != 0 && iSteps < 1000; ++iSteps) {

        // Check for particles with correct anti colour.
        for (int j = iFirst; j < event.size(); j++) {
          if (event[j].isFinal() && !inChain[j] && event[j].acol() == curCol) {
            chain.push_back(j);
            inChain[j] = true;
            curCol = event[j].col();
            break;
          }
        }
        
        // Check for junction with correct colour.
        for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {
          for (int j = 0; j < 3; ++j) {
            if (event.colJunction(iJun,j) == curCol) {
              isJun[isJun.size() -1] = true;
              curCol = 0;
              chain.push_back( -(10 + 10 * iJun + j) );
            }
          }
        }
      }
      chains.push_back(chain);
    }
  }

  // Start from anti-junction and make chains.
  for (int i = 0; i < event.sizeJunction(); ++i) {

    // First check if junction belongs to the right diffractive system.
    int checkCol = event.colJunction(i,0);
    bool wrongSystem = false;
    for (int j = 0; j < iFirst; ++j)
      if (event[j].isFinal() && event[j].acol() == checkCol)
        wrongSystem = true;
    if (wrongSystem)
      continue;

    // Loop over legs of anti junction.
    if (event.kindJunction(i) == 2)
    for (int jCol = 0; jCol < 3; ++jCol) {
      int curCol = event.colJunction(i,jCol);
      vector<int> chain;
      chain.push_back( -(10 + 10 * i + jCol));
      isAntiJun.push_back(true);
      isJun.push_back(false);
      isGluonLoop.push_back(false);
      for (int iSteps = 0; curCol != 0 && iSteps < 1000; ++iSteps) {

        // Check for particles with correct anti colour.
        for (int j = iFirst; j < event.size(); ++j)
        if (event[j].isFinal() && !inChain[j] &&
            event[j].acol() == curCol) {
          chain.push_back(j);
          inChain[j] = true;
          curCol = event[j].col();
          break;
        }
        
        // Check for junction with correct colour.
        for (int iJun = 0; iJun < event.sizeJunction(); ++iJun)
        if (event.kindJunction(iJun) == 1)
        for (int j = 0; j < 3; ++j)
        if (event.colJunction(iJun,j) == curCol) {
          isJun[isJun.size() - 1] = true;
          curCol = 0;
          chain.push_back( -(10 + 10 * iJun + j));
        }
      }
      chains.push_back(chain);
    }
  }

  // Find all gluon loops.
  for (int i = iFirst; i < event.size(); ++i)
  if (event[i].isFinal() && !inChain[i] && event[i].col() != 0) {
    int curCol = event[i].col();
    inChain[i] = true;
    vector<int> chain;
    chain.push_back(i);
    isAntiJun.push_back(false);
    isJun.push_back(false);
    isGluonLoop.push_back(true);
    for (int iSteps = 0; curCol != 0 && iSteps < 1000; ++iSteps) {
      bool foundNext = false;
      for (int j = iFirst; j < event.size(); ++j)
      if (event[j].isFinal() && !inChain[j] && event[j].acol() == curCol) {
        chain.push_back(j);
        inChain[j] = true;
        curCol = event[j].col();
        foundNext = true;
        break;
      }
        
      if (!foundNext)
        break;
    }
    chains.push_back(chain);
  }

  // Form dipoles from chains.
  for (int i = 0; i < int(chains.size()); ++i) {
    if (chains[i].size() == 1) continue;
    int lastCol = -1;
    int firstCol = 0;

    // Start from the first and form the dipoles.
    // Two consececutive dipoles can not share the same colour.
    for (int j = 0; j < int(chains[i].size()); ++j) {
      if (j != int(chains[i].size() - 1)) {

        // Start by picking new colour.
        int newCol = int(rndmPtr->flat() * nReconCols);
        while (newCol == lastCol && !sameNeighbourCol) {
          newCol = int(rndmPtr->flat() * nReconCols);
        }

        // Need to check whether the quark comes from a g->qqbar split.
        // If that is the case, it can not have the same as qbar.
        if (j == 0 && !isAntiJun[i] && !isGluonLoop[i]) {

          int iMother = event[event[ chains[i][j] ].iTopCopy()].mother1();
          if ( event[iMother].idAbs() == 21) {
            vector<int> sisters = event[ chains[i][j] ].sisterList(true);
            // Need to have only one sister and need to be the anti particle.
            if (sisters.size() == 1 && event[ chains[i][j] ].id()
                == - event[ sisters[0] ].id()) {

              // Try to find dipole with sister.
              int colSis = -1;
              for (int k = 0; k < int(dipoles.size()); ++k)
                if (dipoles[k]->iAcol == sisters[0]) {
                  colSis = dipoles[k]->colReconnection;
                  break;
                }
        
              // If the two colours are the same, pick a new.
              while (colSis == newCol && !sameNeighbourCol)
                newCol = int(rndmPtr->flat() * nReconCols);
            }
          }
        }

        // Check if quark end comes from g->qqbar split.
        // If so check that the two quarks get different colours.
        if (j == int(chains[i].size() - 2) && !isJun[i] && !isGluonLoop[i]) {

          int iMother = event[event[chains[i][j + 1]].iTopCopy()].mother1();
          if (event[iMother].idAbs() == 21) {
            vector<int> sisters = event[ chains[i][j + 1] ].sisterList(true);
            // Need to have only one sister and need to be the anti particle.
            if (sisters.size() == 1 && event[ chains[i][j + 1] ].id()
                == - event[ sisters[0] ].id()) {

              // Try to find dipole with sister.
              int colSis = -1;
              for (int k = 0; k < int(dipoles.size()); ++k)
                if (dipoles[k]->iCol == sisters[0]) {
                  colSis = dipoles[k]->colReconnection;
                  break;
                }
        
              // If the two colours are the same, pick a new.
              while ((colSis == newCol || newCol == lastCol)
                     && !sameNeighbourCol)
                newCol = int(rndmPtr->flat() * nReconCols);
            }
          }
        }

        // Special case for junction splitting if multiple gluons
        // between the junctions.
        if ((chains[i][j] > 0 && event[chains[i][j]].status() == 75) ||
            (chains[i][j + 1] > 0 &&
             event[ chains[i][j + 1] ].status() == 75) ) {

          // Find sisters.
          vector<int> sisters;
          if (chains[i][j] > 0 && event[ chains[i][j] ].status() == 75)
            sisters = event[ chains[i][j] ].sisterList();
          else
            sisters = event[ chains[i][j + 1] ].sisterList();

          if (sisters.size() == 3 ) {
        
            // Find colour of sisters.
            int acolSis1 = -1, acolSis2 = -1, acolSis3 = -1;
            int colSis1 = -1, colSis2 = -1, colSis3 = -1;
            for (int k = 0;k < int(dipoles.size()); ++k)  {
              if (dipoles[k]->iAcol == sisters[0])
                acolSis1 = dipoles[k]->colReconnection;

              if (dipoles[k]->iAcol == sisters[1])
                acolSis2 = dipoles[k]->colReconnection;
        
              if (dipoles[k]->iAcol == sisters[2])
                acolSis3 = dipoles[k]->colReconnection;

              if (dipoles[k]->iCol == sisters[0])
                colSis1 = dipoles[k]->colReconnection;

              if (dipoles[k]->iCol == sisters[1])
                colSis2 = dipoles[k]->colReconnection;
        
              if (dipoles[k]->iCol == sisters[2])
                colSis3 = dipoles[k]->colReconnection;
            }

            // If the two colours are the same, pick a new.
            while ((colSis1 == newCol || colSis2 == newCol ||
                   colSis3 == newCol || acolSis1 == newCol ||
                   acolSis2 == newCol || acolSis3 == newCol)
                   && !sameNeighbourCol)
              newCol = int(rndmPtr->flat() * nReconCols);
          }
        }

        // Update stored colours.
        if (j == 0) firstCol = newCol;
        lastCol = newCol;

        // Check if it is anti junction need special dipole.
        if (j == 0 && isAntiJun[i]) {
          int col = event.colJunction( - int(chains[i][j]/10) - 1,
                                       -chains[i][j] % 10);
          dipoles.push_back(new ColourDipole(col, chains[i][j],
            chains[i][j+1], newCol));
          dipoles.back()->isAntiJun = true;
        }

        // Otherwise just make the dipole.
        else dipoles.push_back(new ColourDipole(event[ chains[i][j] ].col(),
          chains[i][j], chains[i][j+1], newCol));

        // If the chain in end a junction mark it.
        if (j == int(chains[i].size() - 2) && isJun[i])
          dipoles.back()->isJun = true;

        // Update the links between dipoles.
        if (j > 0) {
          dipoles[dipoles.size() - 1]->leftDip  = dipoles[dipoles.size() - 2];
          dipoles[dipoles.size() - 2]->rightDip = dipoles[dipoles.size() - 1];
        }
      }

      // If last particle has anti colour it should be possible to connect it
      // to the first particle in the chain. (e.g. gluon loop)
      else
      if (isGluonLoop[i])
      if (event[ chains[i][j] ].col() == event[ chains[i][0] ].acol()) {
        int newCol = int(rndmPtr->flat() * nReconCols);
        while ( (newCol == lastCol || newCol == firstCol)
                && !sameNeighbourCol) {
          newCol = int(rndmPtr->flat() * nReconCols);
        }
        dipoles.push_back(new ColourDipole(event[ chains[i][j] ].col(),
          chains[i][j], chains[i][0], newCol));
        
        // Update links between dipoles.
        dipoles[dipoles.size() - 1]->leftDip = dipoles[dipoles.size() - 2];
        dipoles[dipoles.size() - 2]->rightDip = dipoles[dipoles.size() - 1];
        dipoles[dipoles.size() - chains[i].size()]->leftDip =
          dipoles[dipoles.size() -1];
        dipoles[dipoles.size() - 1]->rightDip =
          dipoles[dipoles.size() - chains[i].size()];
        
      }
    }
  }

  // Setup junction list.
  iColJun.clear();
  iColJun.resize(event.sizeJunction());
  for (int i = 0; i < int(iColJun.size()); ++i) iColJun[i] = vector<int>(3,0);

  // Loop over event and store indices.
  for (int i = iFirst; i < event.size(); ++i)
  if (event[i].isFinal())
  for (int j = 0; j < event.sizeJunction(); ++j)
  for (int jLeg = 0; jLeg < 3; ++jLeg)
  if (event[i].col() == event.colJunction(j,jLeg) ||
      event[i].acol() == event.colJunction(j,jLeg))
    iColJun[j][jLeg] = i;

  // Loop over junction and store indices.
  for (int i = 0;i < event.sizeJunction(); ++i)
  for (int iLeg = 0; iLeg < 3; ++iLeg)
  for (int j = i + 1;j < event.sizeJunction(); ++j)
  for (int jLeg = 0; jLeg < 3; ++jLeg)
  if (event.colJunction(i, iLeg) == event.colJunction(j, jLeg)) {
    iColJun[i][iLeg] = -(10*j + 10 + jLeg);
    iColJun[j][jLeg] = -(10*i + 10 + iLeg);
  }

  // Done.
}

//--------------------------------------------------------------------------

// Calculate the string length of a dipole.

double ColourReconnection::calculateStringLength(ColourDipole * dip,
  vector<ColourDipole *> &dips) {

  if (!dip->isJun && !dip->isAntiJun)  {
    return calculateStringLength(dip->iCol, dip->iAcol);
  }
  else {

    // Start by finding all particles connected to the junction system.
    vector<int> iParticles;
    vector<bool> usedJuns(junctions.size(),false);
    int nJuns = 0;

    if (dip->isJun) {
      if (!findJunctionParticles( -int(dip->iAcol/10) - 1, iParticles,
                                  usedJuns, nJuns, dips)) return 1e9;
    } else
      if (!findJunctionParticles(-int(dip->iCol/10) - 1, iParticles,
                                 usedJuns, nJuns, dips)) return 1e9;

    // If it is a single junction.
    if (int(iParticles.size()) == 3)
      return calculateJunctionLength(iParticles[0], iParticles[1],
                                     iParticles[2]);

    // If it is a junction pair.
    else if (int(iParticles.size()) == 4) {
      return calculateDoubleJunctionLength(iParticles[0], iParticles[1],
                                           iParticles[2], iParticles[3]);
    }
    // If any other number of junction legs return high number.
    else return 1e9;
  }

}
//--------------------------------------------------------------------------

// Update all colours in the event.

void ColourReconnection::updateEvent( Event& event, int iFirst) {

  // Start by making a new copy of particles.
  int oldSize = event.size();
  for (int i = iFirst; i < oldSize;++i)
    if (event[i].isFinal()) event.copy(i, 79);

  // Copy over junctions.
  event.clearJunctions();
  for (int i = 0; i < int(junctions.size()); ++i) {
    for (int j = 0; j < 3; ++j) {
      if ( junctions[i].dipsOrig[j] != 0) {
        junctions[i].col(j, junctions[i].dipsOrig[j]->col);
      }
    }
    event.appendJunction(Junction(junctions[i]));
  }

  // Assign colour according to the real dipoles.
  for (int i = 0; i < int(dipoles.size()); ++i)
    if (dipoles[i]->isReal) {
      if (dipoles[i]->iCol >= 0)
        event[ event[ dipoles[i]->iCol ].daughter1() ].col(dipoles[i]->col);
      else
        event.colJunction(-(dipoles[i]->iCol/10 + 1), -dipoles[i]->iCol % 10,
          dipoles[i]->col);
      if (dipoles[i]->iAcol >= 0)
        event[ event[ dipoles[i]->iAcol ].daughter1() ].acol(dipoles[i]->col);
      else
        event.colJunction(-(dipoles[i]->iAcol/10 + 1), -dipoles[i]->iAcol % 10,
          dipoles[i]->col);
    }
}

//--------------------------------------------------------------------------

// Find all the particles connected in the junction.
// If a single junction, the size of iParticles should be 3.
// For multiple junction structures, the size will increase.

bool ColourReconnection::findJunctionParticles(int iJun,
  vector<int>& iParticles, vector<bool> &usedJuns, int &nJuns,
  vector<ColourDipole*> &dips ) {

  // Mark current junction as used.
  usedJuns[iJun] = true;
  nJuns++;

  // It is not possible to handle junction structures larger than two.
  if (nJuns > 2)
    return false;

  // Find particles connected to the
  if (junctions[iJun].kind() % 2 == 1)
    for (int i = 0; i < 3; ++i)
      iParticles.push_back(junctions[iJun].dips[i]->iCol);
  else
    for (int i = 0; i < 3; ++i)
      iParticles.push_back(junctions[iJun].dips[i]->iAcol);

  // Add dipoles if not already included.
  for (int i = 0; i < 3; ++i) {
    bool addDip = true;
    for (int j = 0; j < int(dips.size()); ++j) {
      if (dips[j] == junctions[iJun].dips[i]) {
        addDip = false;
        break;
      }
    }
    if (addDip) dips.push_back(junctions[iJun].dips[i]);
  }

  // Check whether it connects to any other junctions.
  for (int i = 0; i < int(iParticles.size()); ++i)
  if (iParticles[i] < 0) {
    int iNewJun = - int(iParticles[i] / 10) -1;
    iParticles.erase(iParticles.begin() + i);
    i--;
    if (!usedJuns[iNewJun] && !findJunctionParticles( iNewJun, iParticles,
      usedJuns, nJuns, dips) )
      return false;
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Calculate string length for two indices in the particles record.

double ColourReconnection::calculateStringLength(int i, int j) {
  return stringLength.getStringLength(particles[i].p(), particles[j].p());
}

//--------------------------------------------------------------------------

// Calculate the length of a single junction given the 3 entries in the event.

double ColourReconnection::calculateJunctionLength(int i,
  int j, int k) {

  // Need to be separate indices.
  if ( i == j || i == k || j == k) return 1e9;

  Vec4 p1 = particles[i].p();
  Vec4 p2 = particles[j].p();
  Vec4 p3 = particles[k].p();

  return stringLength.getJuncLength(p1, p2, p3);

}

//--------------------------------------------------------------------------

// Calculate the length of a double junction given the 4 particle entries.
// The first two are expected to be quarks, the second two to be antiquarks.

double ColourReconnection::calculateDoubleJunctionLength( int i, int j,
  int k, int l) {

  // Need to be separate indices.
  if (i == j || i == k || i == l || j == k || j == l || k == l) return 1e9;

  Vec4 p1 = particles[i].p();
  Vec4 p2 = particles[j].p();
  Vec4 p3 = particles[k].p();
  Vec4 p4 = particles[l].p();

  return stringLength.getJuncLength(p1, p2, p3, p4);

}

//--------------------------------------------------------------------------

// Do a single trial emission.

void ColourReconnection::singleReconnection(ColourDipole* dip1,
      ColourDipole* dip2) {

  // Do nothing if it is the same dipole.
  if (dip1 == dip2) return;

  // No colour reconnection if colours do not match.
  if (dip1->colReconnection != dip2->colReconnection) return;

  // If it is not active return
  if (!dip1->isActive || !dip2->isActive) return;

  // Not possible to connect a gluon with itself.
  if (dip1->iCol == dip2->iAcol || dip1->iAcol == dip2->iCol) return;

  // Calculate the difference in lambda.
  double lambdaDiff = getLambdaDiff(dip1, dip2);

  // Insert into trial reconnection if anything is gained.
  if (lambdaDiff > minimumGain) {
    TrialReconnection dipTrial(dip1, dip2, 0, 0, 5, lambdaDiff);
    dipTrials.insert(lower_bound(dipTrials.begin(), dipTrials.end(),
         dipTrial, cmpTrials), dipTrial);
  }
}

//--------------------------------------------------------------------------

// Simple test swap between two dipoles.

void ColourReconnection::swapDipoles(ColourDipole* dip1,
  ColourDipole* dip2, bool back) {

  // Swap the anti colour of the dipoles.
  swap(dip1->iAcol, dip2->iAcol);
  swap(dip1->isJun, dip2->isJun);
  swap(dip1->iAcolLeg, dip2->iAcolLeg);

  // Update the active dipoles. Only change 1 active dipole;
  // this is to avoid problems when switching back.
  if (dip1->iAcol != dip2->iAcol) {
    if (!back) {
      if (dip1->iAcol >= 0)
      for (int i = 0; i < int(particles[dip1->iAcol].activeDips.size()); ++i)
      if (particles[dip1->iAcol].activeDips[i] == dip2) {
        particles[dip1->iAcol].activeDips[i] = dip1;
        swap1 = i;
        break;
      }
      if (dip2->iAcol >= 0)
      for (int i = 0; i < int(particles[dip2->iAcol].activeDips.size()); ++i)
      if (particles[dip2->iAcol].activeDips[i] == dip1) {
        particles[dip2->iAcol].activeDips[i] = dip2;
        swap2 = i;
        break;  
      }
    } else {
      if (dip1->iAcol >= 0) particles[dip1->iAcol].activeDips[swap2] = dip1;
      if (dip2->iAcol >= 0) particles[dip2->iAcol].activeDips[swap1] = dip2;
    }
  }

  // Update list of junctions (only junctions, anti junctions stay the same).
  for (int i = 0; i < int(junctions.size()); ++i)
  if (junctions[i].kind() % 2 == 1)
  for (int iLeg = 0; iLeg < 3; ++iLeg) {
    if (junctions[i].dips[iLeg] == dip1) {
      junctions[i].dips[iLeg] = dip2;
      continue;
    }
    if (junctions[i].dips[iLeg] == dip2) {
      junctions[i].dips[iLeg] = dip1;
      continue;
    }
  }

  // Done.
}

//--------------------------------------------------------------------------

void ColourReconnection::singleJunction(ColourDipole* dip1,
  ColourDipole* dip2) {

   // Do nothing if it is the same dipole.
  if (dip1 == dip2)
    return;

  int iCol1  = dip1->iCol;
  int iCol2  = dip2->iCol;
  int iAcol1 = dip1->iAcol;
  int iAcol2 = dip2->iAcol;

  // Not possible to connet a gluon with itself.
  if (iCol1 == iCol2) return;
  if (iAcol1 == iAcol2) return;

  // Check that all dipoles are active.
  if (!dip1->isActive || !dip2->isActive) return;

  // Do nothing if one of the dipoles is a junction or anti junction.
  if (dip1->isJun || dip1->isAntiJun) return;
  if (dip2->isJun || dip2->isAntiJun) return;

  // Do nothing if it is a pseudo particle that already contains a
  // baryon number inside of it.
  if (int(particles[iCol1].dips.size()) != 1  ||
      int(particles[iAcol1].dips.size()) != 1 ||
      int(particles[iCol2].dips.size()) != 1  ||
      int(particles[iAcol2].dips.size()) != 1)
    return;

  // Only accept 2/9 of the colour configurations.
  if ( (dip1->colReconnection) % 3 !=
        dip2->colReconnection % 3) return;

  if ( (dip1->colReconnection) ==
       dip2->colReconnection) return;

  // Find the colour of the last junction leg.
  int junCol = (3 - (dip1->colReconnection / 3)
    - (dip2->colReconnection / 3) ) * 3
    + (dip1->colReconnection % 3);

  // if other than 9 colours.
  if (nReconCols != 9) {
    while (junCol < 0 || junCol % 3 != dip1->colReconnection % 3 ||
           junCol == dip1->colReconnection || junCol == dip2->colReconnection)
      junCol = int(nReconCols * rndmPtr->flat());
  }

  // Store two new dipoles, these will form the anti-junction.
  ColourDipole* dip3 = dip1;
  ColourDipole* dip4 = dip2;

  double lambdaDiff = getLambdaDiff(dip1, dip2, dip3, dip4, 0);
  if (lambdaDiff > minimumGainJun) {
    TrialReconnection junTrial(dip1, dip2, dip3, dip4, 0, lambdaDiff);
    junTrials.insert(lower_bound(junTrials.begin(), junTrials.end(),
         junTrial, cmpTrials), junTrial);
  }
  // Outer loop
  while (true) {

    // Reset dip4.
    dip4 = dip2;

    // If the colour matches that of the junction.
    if (dip3->colReconnection == junCol)
    while (true) {
      // Check if the new colour matches.
      if (dip4->colReconnection == dip2->colReconnection) {

        // Calculate lambda measure and store new dipole if anything is gained.
        lambdaDiff = getLambdaDiff(dip1, dip2, dip3, dip4, 1);

        if (lambdaDiff > minimumGainJun) {

          TrialReconnection junTrial(dip1, dip2, dip3, dip4, 1, lambdaDiff);
          junTrials.insert(lower_bound(junTrials.begin(), junTrials.end(),
             junTrial, cmpTrials), junTrial);
        }
      }


      // Find the next neighbour.
      if (!findAntiNeighbour(dip4))
        break;

      // Check for gluon loop.
      if (dip4 == dip2 || dip4 == dip1)
        break;

    } // Done with inner loop.

    // Reset dip4.
    dip4 = dip2;

    // If the colour matches that of the other dipole.
    if (dip3->colReconnection == dip1->colReconnection)
    while (true) {
      // Check if the new colour matches.
      if (dip4->colReconnection == junCol) {

        // Calculate lambda measure and store new dipole if anything is gained.
        lambdaDiff = getLambdaDiff(dip1, dip2, dip3, dip4, 2);

        if (lambdaDiff > minimumGainJun) {

          TrialReconnection junTrial(dip1, dip2, dip3, dip4, 2, lambdaDiff);
          junTrials.insert(lower_bound(junTrials.begin(), junTrials.end(),
             junTrial, cmpTrials), junTrial);
        }
      }

      // Find the next neighbour.
      if (!findAntiNeighbour(dip4))
        break;

      // Check for gluon loop.
      if (dip4 == dip2 || dip4 == dip1)
        break;

    } // Done with inner loop.

    // Find the next neighbour.
    if (!findAntiNeighbour(dip3))
      break;

    // Check for gluon loop.
    if (dip3 == dip1 || dip3 == dip2)
      break;
  }

  // Done.

}



//--------------------------------------------------------------------------

void ColourReconnection::singleJunction(ColourDipole* dip1,
  ColourDipole* dip2, ColourDipole* dip3) {

  // Do nothing if one of the dipoles is a junction or anti junction.
  if (dip1->isJun || dip1->isAntiJun) return;
  if (dip2->isJun || dip2->isAntiJun) return;
  if (dip3->isJun || dip3->isAntiJun) return;


  // Check that all dipoles are active.
  if (!dip1->isActive || !dip2->isActive || !dip3->isActive) return;

  // Only allow 0-3-6, 1-4-7 or 2-5-8.
  if ( dip1->colReconnection % 3 != dip2->colReconnection % 3
    || dip1->colReconnection % 3 != dip3->colReconnection % 3) return;

  if ( !(dip1->colReconnection != dip2->colReconnection
      && dip1->colReconnection != dip3->colReconnection
      && dip2->colReconnection != dip3->colReconnection) )
    return;


  if (int(particles[dip1->iCol].dips.size()) != 1  ||
      int(particles[dip1->iAcol].dips.size()) != 1 ||
      int(particles[dip2->iCol].dips.size()) != 1  ||
      int(particles[dip2->iAcol].dips.size()) != 1 ||
      int(particles[dip3->iCol].dips.size()) != 1  ||
      int(particles[dip3->iAcol].dips.size()) != 1 )
    return;

  double lambdaDiff = getLambdaDiff(dip1, dip2, dip3, 0, 3);

  if (lambdaDiff > minimumGainJun) {
    TrialReconnection junTrial(dip1, dip2, dip3, 0, 3, lambdaDiff);
    junTrials.insert(lower_bound(junTrials.begin(), junTrials.end(), junTrial,
                                 cmpTrials), junTrial);
  }

  // Done.
  return;
}

// ------------------------------------------------------------------

void ColourReconnection::makePseudoParticle(ColourDipole* dip , int status,
  bool setupDone) {

  // If it is a normal dipole that needs to be combined.
  if (!dip->isJun && !dip->isAntiJun) {

    // Start by storing variables for easier use.
    int iCol = dip->iCol;
    int iAcol = dip->iAcol;
    int iColLeg = dip->iColLeg;
    int iAcolLeg = dip->iAcolLeg;

    // Make new pseudo particle.
    int iNew = particles.size();
    particles.push_back(particles[iCol]);
    particles[iNew].acol(particles[iCol].acol());
    particles[iNew].col(particles[iAcol].col());
    particles[iNew].mother1(iCol);
    particles[iNew].mother2(iAcol);
    particles[iNew].id(99);
    particles[iNew].status(status);
    particles[iNew].isJun = false;
    particles[iAcol].statusNeg();
    particles[iAcol].daughter1(iNew);
    particles[iCol].statusNeg();
    particles[iCol].daughter1(iNew);
    if (iCol != iAcol)
      particles[iNew].p(particles[iCol].p() + particles[iAcol].p());
    else
      particles[iNew].p(particles[iCol].p());

    // Add all the dipoles from the old pseudo particle.
    // First from particle 1.
    particles[iNew].dips = particles[dip->iCol].dips;
    particles[iNew].colEndIncluded = particles[dip->iCol].colEndIncluded;
    particles[iNew].acolEndIncluded = particles[dip->iCol].acolEndIncluded;

    // Then particle 2.
    if (iCol != iAcol) {
      for (int i = 0; i < int(particles[dip->iAcol].dips.size()); ++i)  {
        if (i != dip->iAcolLeg) {
          // If it is not the same leg, add as separate vector.
          particles[iNew].dips.push_back(particles[dip->iAcol].dips[i]);
          particles[iNew].colEndIncluded.push_back(
            particles[dip->iAcol].colEndIncluded[i]);
          particles[iNew].acolEndIncluded.push_back(
            particles[dip->iAcol].acolEndIncluded[i]);
        } // If it is the same leg, at at the end of the vector.
        else {
          particles[iNew].acolEndIncluded[iColLeg] =
            particles[iAcol].acolEndIncluded[i];
          particles[iNew].dips[iColLeg].pop_back();
          particles[iNew].dips[iColLeg].insert(
            particles[iNew].dips[iColLeg].end(),
            particles[iAcol].dips[i].begin(), particles[iAcol].dips[i].end() );
        }
      }
    }
    if (iCol != iAcol) {
      // Update the dipole legs to the new particle.
      for (int i = 0; i < int(particles[iAcol].activeDips.size()); ++i) {
        if ( particles[iAcol].activeDips[i]->iAcol == iAcol) {
          if (particles[iAcol].activeDips[i]->iAcolLeg < iAcolLeg)
            particles[iAcol].activeDips[i]->iAcolLeg +=
              particles[iCol].dips.size();
          else if (particles[iAcol].activeDips[i]->iAcolLeg == iAcolLeg)
            particles[iAcol].activeDips[i]->iAcolLeg = iColLeg;
          else if (particles[iAcol].activeDips[i]->iAcolLeg > iAcolLeg)
            particles[iAcol].activeDips[i]->iAcolLeg +=
              particles[iCol].dips.size() - 1;
        }
        if (particles[iAcol].activeDips[i]->iCol == iAcol) {
          if (particles[iAcol].activeDips[i]->iColLeg < iAcolLeg)
            particles[iAcol].activeDips[i]->iColLeg +=
              particles[iCol].dips.size();
          else if (particles[iAcol].activeDips[i]->iColLeg == iAcolLeg)
            particles[iAcol].activeDips[i]->iColLeg = iColLeg;
          else if (particles[iAcol].activeDips[i]->iColLeg > iAcolLeg)
            particles[iAcol].activeDips[i]->iColLeg +=
              particles[iCol].dips.size() - 1;
        }
      }
    }

    // Update list of active dipoles.
    particles[iNew].activeDips.clear();
    particles[iNew].activeDips.insert(particles[iNew].activeDips.end(),
      particles[iCol].activeDips.begin(), particles[iCol].activeDips.end());
    if (iCol != iAcol)
    particles[iNew].activeDips.insert(particles[iNew].activeDips.end(),
      particles[iAcol].activeDips.begin(), particles[iAcol].activeDips.end());

    // Remove the now inactive dipole.
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i)
      if (particles[iNew].activeDips[i] == dip) {
        particles[iNew].activeDips.erase(
          particles[iNew].activeDips.begin() + i);
        i--;
      }

    // Update the indices in the active dipoles.
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i) {
      if (particles[iNew].activeDips[i]->iCol == iAcol)
        particles[iNew].activeDips[i]->iCol = iNew;
      if (particles[iNew].activeDips[i]->iCol == iCol)
        particles[iNew].activeDips[i]->iCol = iNew;
      if (particles[iNew].activeDips[i]->iAcol == iAcol)
        particles[iNew].activeDips[i]->iAcol = iNew;
      if (particles[iNew].activeDips[i]->iAcol == iCol)
        particles[iNew].activeDips[i]->iAcol = iNew;
      particles[iNew].activeDips[i]->p1p2
        = mDip(particles[iNew].activeDips[i]);
    }

    // If it is a combination of the same particle,
    // check if any double active dipoles
    if (iCol == iAcol)
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i)
    for (int j = i + 1; j < int(particles[iNew].activeDips.size()); ++j)
    if (particles[iNew].activeDips[i] == particles[iNew].activeDips[j]) {
      particles[iNew].activeDips.erase(particles[iNew].activeDips.begin() + j);
      j--;
    }

    // Add dips changed to used dips.
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i) {
      if (particles[iNew].activeDips[i]->iCol >= 0)
        usedDipoles.push_back(particles[iNew].activeDips[i]);
      else
        for (int j = 0;j < 3; ++j)
          usedDipoles.push_back(junctions[-(particles[iNew].
            activeDips[i]->iCol / 10 + 1)].getColDip(j));

      if (particles[iNew].activeDips[i]->iAcol >= 0)
        usedDipoles.push_back(particles[iNew].activeDips[i]);
      else
        for (int j = 0;j < 3; ++j)
          usedDipoles.push_back(junctions[-(particles[iNew].
            activeDips[i]->iAcol / 10 + 1)].getColDip(j));
    }   

    // mark the internal dipole as not active.
    dip->isActive = false;

    // Done.
    return;
  }

  // If both ends are connected to a junction something went wrong!
  else if (dip->isJun && dip->isAntiJun) {
    return;
  }
  else {

    // Find junction index and first leg to combine.
    int iJun, i0, i1, i2, junLeg0, junLeg1, junLeg2;
    getJunctionIndices(dip, iJun, i0, i1, i2, junLeg0, junLeg1, junLeg2);
    ColourDipole* dip2 = junctions[iJun].dips[junLeg1];
    ColourDipole* dip3 = junctions[iJun].dips[junLeg2];

    // Add new particle.
    int iNew = particles.size();
    particles.push_back(ColourParticle( Particle( 99, status, i0, i1, 0, 0, 0,
      0, particles[i0].p() + particles[i1].p() ) ) );
    particles[iNew].isJun = true;
    particles[iNew].junKind = junctions[iJun].kind();
    if (i0 == i1) particles[iNew].p(particles[i0].p());

    // Update old particles.
    particles[i0].statusNeg();
    particles[i0].daughter1(iNew);
    particles[i1].statusNeg();
    particles[i1].daughter1(iNew);

    // Update list of internal dipoles.
    particles[iNew].dips.clear();
    particles[iNew].dips.insert(particles[iNew].dips.end(),
      particles[i0].dips.begin(),particles[i0].dips.end());
    if (i0 != i1)
      particles[iNew].dips.insert(particles[iNew].dips.end(),
        particles[i1].dips.begin(),particles[i1].dips.end());

    // Update list of whether colour ending is included.
    particles[iNew].colEndIncluded.clear();
    particles[iNew].colEndIncluded.insert(
      particles[iNew].colEndIncluded.end(),
      particles[i0].colEndIncluded.begin(),
      particles[i0].colEndIncluded.end() );
    if (i0 != i1)
      particles[iNew].colEndIncluded.insert(
        particles[iNew].colEndIncluded.end(),
        particles[i1].colEndIncluded.begin(),
        particles[i1].colEndIncluded.end() );

    // Update list of whether anti colour ending is included.
    particles[iNew].acolEndIncluded.clear();
    particles[iNew].acolEndIncluded.insert(
      particles[iNew].acolEndIncluded.end(),
      particles[i0].acolEndIncluded.begin(),
      particles[i0].acolEndIncluded.end() );
    if (i0 != i1)
      particles[iNew].acolEndIncluded.insert(
        particles[iNew].acolEndIncluded.end(),
        particles[i1].acolEndIncluded.begin(),
        particles[i1].acolEndIncluded.end() );

    // Third particle just need to add one to list of dipoles.
    if (dip->isJun && i2 >= 0 && i2 != i0 && i2 != i1) {
      particles[iNew].dips.push_back(particles[i2].dips[dip3->iColLeg]);
      particles[iNew].dips.back().erase(particles[iNew].dips.back().begin(),
        particles[iNew].dips.back().end() - 1);

    }
    if (dip->isAntiJun && i2 >= 0 && i2 != i0 && i2 != i1) {
      particles[iNew].dips.push_back(particles[i2].dips[dip3->iAcolLeg]);
      particles[iNew].dips.back().erase(
        particles[iNew].dips.back().begin() + 1,
        particles[iNew].dips.back().end() );
    }

    // Add endings for the third particle.
    if (i2 != i0 && i2 != i1) {
      particles[iNew].acolEndIncluded.push_back(false);
      particles[iNew].colEndIncluded.push_back(false);
    }

    // Special case if it is J-J connection.
    if (i2 < 0) {
      particles[iNew].dips.push_back(vector<ColourDipole *>());

      // Find the real dipole to add to dipole list.
      for (int i = 0; i < int(dipoles.size()); ++i)
        if (dipoles[i]->isReal && dipoles[i]->iCol == dip3->iCol &&
            dipoles[i]->iAcol == dip3->iAcol)
          particles[iNew].dips.back().push_back(dipoles[i]);

      // Change ending.
      particles[iNew].acolEndIncluded.back() = true;
      particles[iNew].colEndIncluded.back()  = true;
    }

    // The endings need to reflect the new junction structure.
    if (dip->isJun)
    for (int i = 0; i < int(particles[iNew].acolEndIncluded.size()); ++i)
      particles[iNew].acolEndIncluded[i] = true;
    else
    for (int i = 0; i < int(particles[iNew].colEndIncluded.size()); ++i)
      particles[iNew].colEndIncluded[i] = true;

    // Update active dipoles, first junction case.
    // Set the now internal dipoles as inactive.
    dip->isActive = false;
    dip2->isActive = false;
    dip3->isActive = true;

    // Update the dipole legs to the new particle.
    // Only need to do it for the iAcol particle,
    // since nothing changes for the iCol particle.
    if (i0 != i1)
    for (int i = 0; i < int(particles[i1].activeDips.size()); ++i) {
      if (particles[i1].activeDips[i]->iAcol == i1)
        particles[i1].activeDips[i]->iAcolLeg += particles[i0].dips.size();
      if (particles[i1].activeDips[i]->iCol == i1)
        particles[i1].activeDips[i]->iColLeg += particles[i0].dips.size();
    }

    // Update list of active dipoles.
    particles[iNew].activeDips.clear();
    particles[iNew].activeDips.insert(particles[iNew].activeDips.end(),
      particles[i0].activeDips.begin(), particles[i0].activeDips.end());
    if (i0 != i1)
      particles[iNew].activeDips.insert(particles[iNew].activeDips.end(),
        particles[i1].activeDips.begin(), particles[i1].activeDips.end());
    if (i2 != i0 && i2 != i1)
      particles[iNew].activeDips.push_back(dip3);

    // Remove the now inactive dipoles.
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i) {
      if (particles[iNew].activeDips[i] == dip) {
        particles[iNew].activeDips.erase(
          particles[iNew].activeDips.begin() + i);
        i--;
        continue;
      }
      if (particles[iNew].activeDips[i] == dip2) {
        particles[iNew].activeDips.erase(
          particles[iNew].activeDips.begin() + i);
        i--;
        continue;
      }
    }

    // Update the indices in the active dipoles.
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i) {
      if (particles[iNew].activeDips[i]->iCol == i1)
        particles[iNew].activeDips[i]->iCol = iNew;
      if (particles[iNew].activeDips[i]->iCol == i0)
        particles[iNew].activeDips[i]->iCol = iNew;
      if (particles[iNew].activeDips[i]->iAcol == i1)
        particles[iNew].activeDips[i]->iAcol = iNew;
      if (particles[iNew].activeDips[i]->iAcol == i0)
        particles[iNew].activeDips[i]->iAcol = iNew;
      particles[iNew].activeDips[i]->p1p2
        = mDip(particles[iNew].activeDips[i]);
    }

    // The third dip is no longer connected to a junction.
    if (dip->isJun) {
      dip3->isJun = false;
      dip3->iAcol = iNew;
      if (i2 != i0 && i2 != i1)
        dip3->iAcolLeg = particles[iNew].dips.size() - 1;
    }
    else  {
      dip3->isAntiJun = false;
      dip3->iCol = iNew;
      if (i2 != i0 && i2 != i1)
        dip3->iColLeg = particles[iNew].dips.size() - 1;
    }

    // Add dips changed to used dips.
    for (int i = 0; i < int(particles[iNew].activeDips.size()); ++i) {
      bool added = false;
      for (int j = 0;j < int(usedDipoles.size()); ++j)
        if (particles[iNew].activeDips[i] == usedDipoles[j]) {
          added = true;
          break;
        }
      if (!added) usedDipoles.push_back(particles[iNew].activeDips[i]);
    }
    usedDipoles.push_back(dip);
    usedDipoles.push_back(dip2);
    usedDipoles.push_back(dip3);


    // Possible for the new dip to have a low m0.
    if (setupDone && mDip(dip3) < m0)
      makePseudoParticle(dip3, status, true);
  }

  // Done.

}

// ------------------------------------------------------------------

// Help function to sort dipoles in right order.

bool sortFunc(ColourDipole* a, ColourDipole* b) {
    return (a->p1p2 < b->p1p2);
}

// ------------------------------------------------------------------

// Form all pseudoparticles below m0.

void ColourReconnection::makeAllPseudoParticles( Event & event, int iFirst) {

  // Make junctions.
  for (int i = 0; i < event.sizeJunction(); ++i)
    junctions.push_back(event.getJunction(i));

  // Make new copy of all the dipoles.
  int oldSize = int(dipoles.size());
  for (int i = 0; i < oldSize; ++i) {
    dipoles.push_back(new ColourDipole(*dipoles[i]));
    dipoles[i + oldSize]->iColLeg = 0;
    dipoles[i + oldSize]->iAcolLeg = 0;
    dipoles[i]->iColLeg = 0;
    dipoles[i]->iAcolLeg = 0;
    dipoles[i]->isActive = false;
    dipoles[i]->isReal = true;
    dipoles[i + oldSize]->isReal = false;

    // Store original dipoles connected to junctions.
    if (dipoles[i]->iCol < 0) {
      junctions[-(dipoles[i]->iCol / 10 + 1)].dipsOrig[(-dipoles[i]->iCol)
        % 10] = dipoles[i];
    }
    if (dipoles[i]->iAcol < 0) {
      junctions[-(dipoles[i]->iAcol / 10 + 1)].dipsOrig[-(dipoles[i]->iAcol
        % 10)] = dipoles[i];
    }
  }

  // Set up the coldDips and acolDips.
  for (int i = 0; i < oldSize; ++i) {
    if (dipoles[i]->leftDip != 0)
    for (int j = 0; j < oldSize; ++j)
    if (dipoles[i]->leftDip == dipoles[j]) {
      dipoles[i + oldSize]->colDips.push_back(dipoles[j + oldSize]);
      break;
    }

    if (dipoles[i]->rightDip != 0)
    for (int j = 0; j < oldSize; ++j)
    if (dipoles[i]->rightDip == dipoles[j]) {
      dipoles[i + oldSize]->acolDips.push_back(dipoles[j + oldSize]);
      break;
    }
  }

  // Start by copying event record to make pseudoparticles.
  // The pseudoparticles also need to gain
  for (int i = iFirst; i < event.size(); ++i)
  if (event[i].isFinal()) {
    particles.push_back(ColourParticle(event[i]));
    particles.back().dips.resize(1,vector<ColourDipole *>());

    // Set up dipoles.
    for (int j = 0; j < int(dipoles.size()); ++j) {
      if (dipoles[j]->iCol == i) {
        if (dipoles[j]->isActive) {
          dipoles[j]->iCol = particles.size() - 1;
          particles.back().activeDips.push_back(dipoles[j]);
        }
        else particles.back().dips[0].push_back(dipoles[j]);
      }

      if (dipoles[j]->iAcol == i) {
        if (dipoles[j]->isActive) {
          dipoles[j]->iAcol = particles.size() - 1;
          particles.back().activeDips.push_back(dipoles[j]);
        }
        else particles.back().dips[0].insert(particles.back().dips[0].begin(),
          dipoles[j]);
      }
    }

    // Tell whether dipoles are connected to other dipoles.
    if (event[i].isQuark() && event[i].id() > 0)
      particles.back().colEndIncluded.push_back(true);
    else particles.back().colEndIncluded.push_back(false);

    if (event[i].isQuark() && event[i].id() < 0)
      particles.back().acolEndIncluded.push_back(true);
    else particles.back().acolEndIncluded.push_back(false);
  }

  // Inserting a copy of the event record, but now with full
  // pseudo particle setup.
  // This is mainly to avoid having to distinguish between combining
  // original particles and pseudoparticles.

  // Set right dipole connections in junctions.
  for (int i = 0; i < int(dipoles.size()); ++i) {
    if (dipoles[i]->iCol < 0) {
      int j = (- dipoles[i]->iCol / 10) - 1;
      int jLeg = - dipoles[i]->iCol % 10;
      junctions[j].setColDip(jLeg, dipoles[i]);
    }
    if (dipoles[i]->iAcol < 0) {
      int j = (- dipoles[i]->iAcol / 10) - 1;
      int jLeg = - dipoles[i]->iAcol % 10;
      junctions[j].setColDip(jLeg, dipoles[i]);
    }
  }

  // Make sure all dipoles masses are set correctly.
  for (int i = 0; i < int(dipoles.size()); ++i) {
    if (dipoles[i]->isActive)
      dipoles[i]->p1p2 = mDip(dipoles[i]);
    else
      dipoles[i]->p1p2 = 1e9;
  }

  // Keep making pseudo particles until they are above the threshold.
  while (true) {
    sort(dipoles.begin(), dipoles.end(), sortFunc);
    bool finished = true;
    for (int i = 0; i < int(dipoles.size()); ++i) {
      if (!dipoles[i]->isActive) continue;
      if (dipoles[i]->p1p2 < m0) {
        makePseudoParticle( dipoles[i], 110);
        finished = false;
        break;
      }
      else break;
    }
    if (finished) break;
  }

  // Sort the dipoles.
  sort(dipoles.begin(), dipoles.end(), sortFunc);

  // Done.
  return;

}

// ------------------------------------------------------------------

// Print statements if something is wrong in dipole setup.
// Does not have a return statement -- DEBUG PURPOSE ONLY --.

void ColourReconnection::checkRealDipoles(Event& event, int iFirst) {
  vector<int> dipConnections(event.size(),0);
  for (int i = 0;i < int(dipoles.size()); ++i)
    if (dipoles[i]->isReal) {
      if (dipoles[i]->iCol >= 0)
        dipConnections[dipoles[i]->iCol]++;
      if (dipoles[i]->iAcol >= 0)
        dipConnections[dipoles[i]->iAcol]++;
    }
  bool working = true;
  for (int i = iFirst ;i < event.size(); ++i) {
    if (event[i].isFinal()) {
      if (event[i].isQuark() && dipConnections[i] != 1) {
        cout << "quark " << i << " is wrong!!" << endl;
        working = false;
      }
      else if (event[i].idAbs() == 21 && dipConnections[i] != 2) {
        cout << "gluon " << i << " is wrong!!" << endl;
        working = false;
      }
    }
  }
  if (!working) {
    infoPtr->errorMsg("Error in ColourReconnection::checkRealDipoles:"
      "Real dipoles not setup properply");

  }

}
// ------------------------------------------------------------------

// Print statements if something is wrong in dipole setup.
// Does not have a return statement -- DEBUG PURPOSE ONLY --.

void ColourReconnection::checkDipoles() {

  for (int i = 0;i < int(dipoles.size()); ++i) {
    if (dipoles[i] == 0) { cout << "dipole empty" << endl;}
    if (dipoles[i]->isActive) {
      if (dipoles[i]->iCol >= 0) {
        bool foundMyself = false;
        for (int j = 0; j < int(particles[ dipoles[i]->iCol ].
          activeDips.size()); ++j) {
          if (!particles[dipoles[i]->iCol].activeDips[j]->isActive) {
            infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
              "Found inactive dipole, where only active was expected");
          }
          if (particles[dipoles[i]->iCol].activeDips[j] == dipoles[i])
            foundMyself = true;
        }

        if (!foundMyself) {
          infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
            "Linking between active dipoles and particles is wrong");
        }
        if (dipoles[i]->iColLeg
          >= int(particles[dipoles[i]->iCol].dips.size())) {
          infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
            "Original dipoles not stored correct");
        }

        // Check that linking to old dipoles work.
        if (dipoles[i]->col !=
           particles[dipoles[i]->iCol].dips[dipoles[i]->iColLeg].back()->col) {
           infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
            "Original dipoles do not match in");
        }
      }

      if (dipoles[i]->iAcol >= 0) {
        bool foundMyself = false;
        for (int j = 0;j < int(particles[ dipoles[i]->iAcol ].
          activeDips.size()); ++j) {
        
          if (!particles[dipoles[i]->iAcol].activeDips[j]->isActive) {
            infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
              "Found inactive dipole, where only active was expected");
          }
           if (particles[dipoles[i]->iAcol].activeDips[j] == dipoles[i])
            foundMyself = true;
        }

        if (!foundMyself) {
           infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
            "Linking between active dipoles and particles is wrong");
        }
        if (dipoles[i]->iAcolLeg >= int(particles[dipoles[i]->iAcol].
          dips.size() )) {
          infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
            "Original dipoles not stored correct");
        }

        // Check that linking to old dipoles work
        if (dipoles[i]->col != particles[dipoles[i]->iAcol].
            dips[dipoles[i]->iAcolLeg].front()->col) {
           infoPtr->errorMsg("Error in ColourReconnection::checkDipoles:"
            "Original dipoles do not match in");
        }
      }
    }
  }
}

// ------------------------------------------------------------------

// Print all the chains.

void ColourReconnection::listAllChains() {

  cout << "  ----- PRINTING CHAINS -----  " << dipoles.size() << endl;
  for (int i = 0; i < int(dipoles.size()); ++i)
    dipoles[i]->printed = false;

  for (int i = 0;i < int(dipoles.size()); ++i)
    if (!dipoles[i]->printed)
      listChain(dipoles[i]);
  cout << "  ----- PRINTED CHAINS -----  " << endl;

}

// ------------------------------------------------------------------

// Print the chain containing the dipole.

void ColourReconnection::listChain(ColourDipole *dip) {

  // Make sure not an empty pointer.
  if (dip == 0) return;

  // If chain is not active, just print it.
  if (!dip->isActive) {
    return;
  }

  ColourDipole * colDip = dip;

  // Try to reach one end of the chain.
  while (particles[colDip->iCol].dips.size() == 1 && findColNeighbour(colDip))
    if (dip == colDip)
      break;

  ColourDipole * endDip = colDip;
  do {
    cout << colDip->iCol << " (" << colDip->p1p2 << ", " << colDip->col
         << ") (" << colDip->isActive << ") ";
    colDip->printed = true;
  }
  // Start the printing.
  while (particles[colDip->iAcol].dips.size() == 1 && findAntiNeighbour(colDip)
         && colDip != endDip);

  // Print the last part.
  cout << colDip->iAcol<< endl;

  // Done.
}

// ------------------------------------------------------------------

// Return relevant indices for the junction.

bool ColourReconnection::getJunctionIndices(ColourDipole * dip, int &iJun,
  int &i0, int &i1, int &i2, int &junLeg0, int &junLeg1, int &junLeg2) {

  // Find junction index and first leg to combine.
  int indxJun = dip->iCol;
  if (dip->iAcol < 0)
      indxJun = dip->iAcol;
  iJun = (- indxJun / 10) - 1;
  junLeg0 = -(indxJun % 10);
  junLeg1 = 1;
  junLeg2 = 2;
  if (junLeg0 == 1) junLeg1 = 0;
  else if (junLeg0 == 2) junLeg2 = 0;

  if (dip->iCol < 0) {
    i0 = dip->iAcol;
    i1 = junctions[iJun].dips[junLeg1]->iAcol;
    i2 = junctions[iJun].dips[junLeg2]->iAcol;
  }
  else {
    i0 = dip->iCol;
    i1 = junctions[iJun].dips[junLeg1]->iCol;
    i2 = junctions[iJun].dips[junLeg2]->iCol;
  }

  // It is not possible to form a pseudoparticle if only a single particle is
  // connected to the junction.
  if (i1 < 0 && i2 < 0) return false;

  // Check which two particle should form the pseudoparticle.
  double m1 = 1e9, m2 = 1e9;
  if (i1 >= 0)
    m1 = m(particles[i0].p(),particles[i1].p());
  if (i2 >= 0)
    m2 = m(particles[i0].p(),particles[i2].p());

  if (m1 > m2) {
    swap(i1,i2);
    swap(junLeg1,junLeg2);
  }
  // Force switch if i0 == i2
  if (i0 == i2) {
    swap(i1,i2);
    swap(junLeg1,junLeg2);
  }

  return true;
}

// ------------------------------------------------------------------

// Calculate the invariant mass of a dipole.

double ColourReconnection::mDip(ColourDipole* dip) {

  // If double junction no invariant mass is given.
  if (dip->isJun && dip->isAntiJun) return 1e9;
  // If it has a single junction end.
  else if (dip->isJun || dip->isAntiJun) {
    int iJun, i0, i1, i2, junLeg0, junLeg1, junLeg2;
    getJunctionIndices(dip, iJun, i0, i1, i2, junLeg0, junLeg1, junLeg2);
    if (i0 == i1)
      return particles[i0].m();
    if (i1 < 0)
      return 1e9;
    return m(particles[i0].p(),particles[i1].p());
  } // No junction ends.
  else {
    if (dip->iCol == dip->iAcol)
      return particles[dip->iCol].m();
    else
      return m(particles[dip->iCol].p(),particles[dip->iAcol].p());
  }

}

// ------------------------------------------------------------------

// Print dipoles, intended for debuggning purposes.

void ColourReconnection::listDipoles(bool onlyActive, bool onlyReal) {

  cout << " --- listing dipoles ---" << endl;
  for (int i = 0; i < int(dipoles.size()); ++i) {
    if (onlyActive && !dipoles[i]->isActive)
      continue;
    if (onlyReal && !dipoles[i]->isReal)
      continue;
    dipoles[i]->print();
  }
  cout << " --- finished listing ---" << endl;

}

// ------------------------------------------------------------------

// Print particles, intended for debugging purposes.

void ColourReconnection::listParticles() {

  for (int i = 0; i < int(particles.size()); ++i) {
    const ColourParticle& pt = particles[i];

    // Basic line for a particle, always printed.
    cout << setw(6) << i << setw(10) << pt.id() << "   " << left
       << setw(18) << pt.nameWithStatus(18) << right << setw(4)
       << pt.status() << setw(6) << pt.mother1() << setw(6)
       << pt.mother2() << setw(6) << pt.daughter1() << setw(6)
       << pt.daughter2() << setw(6) << pt.col() << setw(6) << pt.acol()
       << setprecision(3)
       << setw(11) << pt.px() << setw(11) << pt.py() << setw(11)
       << pt.pz() << setw(11) << pt.e() << setw(11) << pt.m();
    for (int j = 0;j < int(pt.activeDips.size());++j)
      cout << setw(10) << pt.activeDips[j];
    cout << "\n";
  }

}

// ------------------------------------------------------------------

// Print junctions, intended for debugging purposes.

void ColourReconnection::listJunctions() {

  cout << " --- listing junctions ---" << endl;
  for (int i = 0; i < int(junctions.size()); ++i)
    junctions[i].print();
  cout << " --- finished listing ---" << endl;

}

// ------------------------------------------------------------------

// Allow colour reconnections by mergings of MPI collision subsystems.
// iRec is system that may be reconnected, by moving its gluons to iSys,
// where minimal pT (or equivalently Lambda) is used to pick location.
// Therefore all dipoles in iSys have to be found, and all gluons in iRec.
// Matching q-qbar pairs are treated by analogy with gluons.
// Note: owing to rescatterings some outgoing partons must be skipped.

bool ColourReconnection::reconnectMPIs( Event&  event, int oldSize) {

  // References to beams to simplify indexing.
  BeamParticle& beamA = *beamAPtr;
  BeamParticle& beamB = *beamBPtr;

  // Prepare record of which systems should be merged onto another.
  // The iSys system must have colour in final state to attach to it.
  nSys = partonSystemsPtr->sizeSys();
  vector<int>  iMerge(nSys);
  vector<bool> hasColour(nSys);
  for (int iSys = 0; iSys < nSys; ++iSys) {
    iMerge[iSys] = iSys;
    bool hasCol = false;
    for (int iMem = 0; iMem < partonSystemsPtr->sizeOut(iSys); ++iMem) {
      int iNow = partonSystemsPtr->getOut( iSys, iMem);
      if (event[iNow].isFinal() && (event[iNow].col() > 0
        || event[iNow].acol() > 0) ) {
        hasCol = true;
        break;
      }
    }
    hasColour[iSys] = hasCol;
  }

  // Loop over systems to decide which should be reconnected.
  for (int iRec = nSys - 1; iRec > 0; --iRec) {

    // Determine reconnection strength from pT scale of system.
    double pT2Rec  = pow2( partonSystemsPtr->getPTHat(iRec) );
    double probRec = pT20Rec / (pT20Rec + pT2Rec);

    // Loop over other systems iSys at higher pT scale and
    // decide whether to reconnect the iRec gluons onto one of them.
    for (int iSys = iRec - 1; iSys >= 0; --iSys)
    if (hasColour[iSys] && probRec > rndmPtr->flat()) {

      // The iRec system and all merged with it to be merged with iSys.
      iMerge[iRec] = iSys;
      for (int iRec2 = iRec + 1; iRec2 < nSys; ++iRec2)
      if (iMerge[iRec2] == iRec) iMerge[iRec2] = iSys;

      // Once a system has been merged do not test it anymore.
      break;
    }
  }

  // Loop over systems. Check whether other systems to be merged with it.
  for (int iSys = 0; iSys < nSys; ++iSys) {
    int nMerge = 0;
    for (int iRec = iSys + 1; iRec < nSys; ++iRec)
    if (iMerge[iRec] == iSys) ++nMerge;
    if (nMerge == 0) continue;

    // Incoming partons not counted if rescattered.
    int  iInASys = partonSystemsPtr->getInA(iSys);
    bool hasInA  = (beamA[iSys].isFromBeam());
    int  iInBSys = partonSystemsPtr->getInB(iSys);
    bool hasInB  = (beamB[iSys].isFromBeam());

    // Begin find dipoles in iSys system.
    vector<BeamDipole> bmdipoles;
    int sizeOut = partonSystemsPtr->sizeOut(iSys);
    for (int iMem = 0; iMem < sizeOut; ++iMem) {

      // Find colour dipoles to beam remnant.
      int iNow = partonSystemsPtr->getOut( iSys, iMem);
      if (!event[iNow].isFinal()) continue;
      int col = event[iNow].col();
      if (col > 0) {
        if      (hasInA && event[iInASys].col() == col)
          bmdipoles.push_back( BeamDipole( col, iNow, iInASys ) );
        else if (hasInB && event[iInBSys].col() == col)
          bmdipoles.push_back( BeamDipole( col, iNow, iInBSys ) );

        // Find colour dipole between final-state partons.
        else for (int iMem2 = 0; iMem2 < sizeOut; ++iMem2)
        if (iMem2 != iMem) {
          int iNow2 = partonSystemsPtr->getOut( iSys, iMem2);
          if (!event[iNow2].isFinal()) continue;
          if (event[iNow2].acol() == col) {
            bmdipoles.push_back( BeamDipole( col, iNow, iNow2) );
            break;
          }
        }
      }

      // Find anticolour dipoles to beam remnant.
      int acol = event[iNow].acol();
      if (acol > 0) {
        if      (hasInA && event[iInASys].acol() == acol)
          bmdipoles.push_back( BeamDipole( acol, iInASys, iNow ) );
        else if (hasInB && event[iInBSys].acol() == acol)
          bmdipoles.push_back( BeamDipole( acol, iInBSys, iNow ) );
      }
    }

    // Skip mergings if no dipoles found.
    if (bmdipoles.size() == 0) continue;

    // Find dipole sizes.
    for (int iDip = 0; iDip < int(bmdipoles.size()); ++iDip)
      bmdipoles[iDip].p1p2 = event[bmdipoles[iDip].iCol].p()
                           * event[bmdipoles[iDip].iAcol].p();

    // Loop over systems iRec to be merged with iSys.
    for (int iRec = iSys + 1; iRec < nSys; ++iRec) {
      if (iMerge[iRec] != iSys) continue;

      // Information on iRec. Vectors for gluons and anything else.
      int sizeRec = partonSystemsPtr->sizeOut(iRec);
      int iInARec = partonSystemsPtr->getInA(iRec);
      int iInBRec = partonSystemsPtr->getInB(iRec);
      int nGluRec = 0;
      vector<int>    iGluRec;
      vector<double> pT2GluRec;
      int nAnyRec = 0;
      vector<int>    iAnyRec;
      vector<bool>   freeAnyRec;

      // Copy of gluon positions in descending order.
      for (int iMem = 0; iMem < sizeRec; ++iMem) {
        int iNow = partonSystemsPtr->getOut( iRec, iMem);
        if (!event[iNow].isFinal()) continue;
        if (event[iNow].isGluon()) {
          ++nGluRec;
          iGluRec.push_back( iNow );
          pT2GluRec.push_back( event[iNow].pT2() );
          for (int i = nGluRec - 1; i > 1; --i) {
            if (pT2GluRec[i - 1] > pT2GluRec[i]) break;
            swap(   iGluRec[i - 1],   iGluRec[i] );
            swap( pT2GluRec[i - 1], pT2GluRec[i] );
          }
        // Copy of anything else, mainly quarks, in no particular order.
        } else {
          ++nAnyRec;
          iAnyRec.push_back( iNow );
          freeAnyRec.push_back( true );
        }
      }

      // For each gluon in iRec now find the dipole that gives the smallest
      // (pGlu * pI) (pGlu * pJ) / (pI * pJ), i.e. minimal pT (and Lambda).
      for (int iGRec = 0; iGRec < nGluRec; ++iGRec) {
        int    iGlu      = iGluRec[iGRec];
        Vec4   pGlu      = event[iGlu].p();
        int    iDipMin   = 0;
        double pT2DipMin = sCM;
        for (int iDip = 0; iDip < int(bmdipoles.size()); ++iDip) {
          double pT2Dip = (pGlu * event[bmdipoles[iDip].iCol].p())
            * (pGlu * event[bmdipoles[iDip].iAcol].p()) / bmdipoles[iDip].p1p2;
          if (pT2Dip < pT2DipMin) {
            iDipMin   = iDip;
            pT2DipMin = pT2Dip;
          }
        }

        // Attach the gluon to the dipole, i.e. split the dipole in two.
        int colGlu   = event[iGlu].col();
        int acolGlu  = event[iGlu].acol();
        int colDip   = bmdipoles[iDipMin].col;
        int iColDip  = bmdipoles[iDipMin].iCol;
        int iAcolDip = bmdipoles[iDipMin].iAcol;
        event[iGlu].acol( colDip );
        if (event[iAcolDip].acol() == colDip)
             event[iAcolDip].acol( colGlu );
        else event[iAcolDip].col(  colGlu );
        bmdipoles[iDipMin].iAcol = iGlu;
        bmdipoles[iDipMin].p1p2 = event[iColDip].p() * pGlu;
        bmdipoles.push_back( BeamDipole( colGlu, iGlu, iAcolDip ) );
        bmdipoles.back().p1p2 = pGlu * event[iAcolDip].p();

        // Remove gluon from old system: reconnect colours.
        for (int i = oldSize; i < event.size(); ++i)
        if (i != iGlu && i != iAcolDip) {
          if (event[i].isFinal()) {
            if (event[i].acol() == colGlu) event[i].acol( acolGlu );
          } else {
              if (event[i].col()  == colGlu) event[i].col( acolGlu );
          }
        }

        // Update any junction legs that match reconnected dipole.
        for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {

          // Only junctions need to be updated, not antijunctions.
          if (event.kindJunction(iJun) % 2 == 0) continue;
          for (int leg = 0; leg < 3; ++leg) {
            int col = event.colJunction(iJun, leg);
            if (col == colDip)
              event.colJunction(iJun, leg, colGlu);
          }
        }

      }

      // See if any matching quark-antiquark pairs among the rest.
      for (int iQRec = 0; iQRec < nAnyRec; ++iQRec) {
        int iQ  = iAnyRec[iQRec];
        int idQ = event[iQ].id();
        if (freeAnyRec[iQRec] && idQ > 0 && idQ < 6)
        for (int iQbarRec = 0; iQbarRec < nAnyRec; ++iQbarRec) {
          int iQbar  = iAnyRec[iQbarRec];
          if (freeAnyRec[iQbarRec] && event[iQbar].id() == -idQ) {

            // Check that these can be traced back to same gluon splitting.
            // For now also avoid qqbar pairs produced in rescatterings.??
            int iTopQ    = event[iQ].iTopCopyId();
            int iTopQbar = event[iQbar].iTopCopyId();
            int iMother  = event[iTopQ].mother1();
            if (event[iTopQbar].mother1() == iMother
              && event[iMother].isGluon() && event[iMother].status() != -34
              && event[iMother + 1].status() != -34 ) {

              // Now find the dipole that gives the smallest
              // ((pQ + pQbar) * pI) ((pQ + pQbar) * pJ) / (pI * pJ).
              Vec4   pGlu      = event[iQ].p() + event[iQbar].p();
              int    iDipMin   = 0;
              double pT2DipMin = sCM;
              for (int iDip = 0; iDip < int(bmdipoles.size()); ++iDip) {
                double pT2Dip = (pGlu * event[bmdipoles[iDip].iCol].p())
                  * (pGlu * event[bmdipoles[iDip].iAcol].p())
                  / bmdipoles[iDip].p1p2;
                if (pT2Dip < pT2DipMin) {
                  iDipMin   = iDip;
                  pT2DipMin = pT2Dip;
                }
              }

              // Attach the q-qbar pair to the dipole, i.e. split the dipole.
              int colGlu   = event[iQ].col();
              int acolGlu  = event[iQbar].acol();
              int colDip   = bmdipoles[iDipMin].col;
              int iColDip  = bmdipoles[iDipMin].iCol;
              int iAcolDip = bmdipoles[iDipMin].iAcol;
              event[iQbar].acol( colDip );
              if (event[iAcolDip].acol() == colDip)
                   event[iAcolDip].acol( colGlu );
              else event[iAcolDip].col(  colGlu );
              bmdipoles[iDipMin].iAcol = iQbar;
              bmdipoles[iDipMin].p1p2 = event[iColDip].p() * event[iQbar].p();
              bmdipoles.push_back( BeamDipole( colGlu, iQ, iAcolDip ) );
              bmdipoles.back().p1p2 = event[iQ].p() * event[iAcolDip].p();

              // Remove q-qbar pair from old system: reconnect colours.
              freeAnyRec[iQRec]    = false;
              freeAnyRec[iQbarRec] = false;
              for (int i = oldSize; i < event.size(); ++i)
              if (i != iQRec && i != iQbarRec && i != iColDip
                && i != iAcolDip) {
                if (event[i].isFinal()) {
                  if (event[i].acol() == colGlu) event[i].acol( acolGlu );
                } else {
                    if (event[i].col()  == colGlu) event[i].col( acolGlu );
                }
              }

              // Update any junction legs that match reconnected dipole.
              for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {

                // Only junctions need to be updated, not antijunctions.
                if (event.kindJunction(iJun) % 2 == 0) continue;
                for (int leg = 0; leg < 3; ++leg) {
                  int col = event.colJunction(iJun, leg);
                  if (col == colDip)
                    event.colJunction(iJun, leg, colGlu);
                }
              }

              // Done with processing of q-qbar pairs.
            }
          }
        }
      }

      // If only two beam gluons left of system, set their colour = anticolour.
      // Used by BeamParticle::remnantColours to skip irrelevant gluons.
      if ( event[iInARec].isGluon() && !event[iInARec].isRescatteredIncoming()
        && event[iInBRec].isGluon() && !event[iInBRec].isRescatteredIncoming()
        && event[iInARec].col() == event[iInBRec].acol()
        && event[iInARec].acol() == event[iInBRec].col() ) {
          event[iInARec].acol( event[iInARec].col() );
          event[iInBRec].acol( event[iInBRec].col() );
      }

    // End of loops over iRec and iSys systems.
    }
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Find the neighbour to the anticolour side. Return false if the dipole
// is connected to a junction or the new particle has a junction inside of it.

bool ColourReconnection::findAntiNeighbour(ColourDipole*& dip) {
  // If only one active dipole, it has to be an antiquark.
  if (int(particles[dip->iAcol].activeDips.size())  == 1)
    return false;

  // Has to have to active dipoles, otherwise something went wrong.
  if (int(particles[dip->iAcol].activeDips.size())  != 2) {
    infoPtr->errorMsg("Warning in ColourReconnection::findAntiNeighbour: "
                      "Wrong number of active dipoles");
    return false;
  }

  // Otherwise find new dipole.
  if (dip == particles[dip->iAcol].activeDips[0])
    dip = particles[dip->iAcol].activeDips[1];
  else dip = particles[dip->iAcol].activeDips[0];

  // Do not allow the new dipole to be connected to an antijunction.
  if (dip->isAntiJun || dip->isJun)
    return false;

  // Do not allow new dipole to have a pseudoparticle with
  // a baryon number inside.
  if (int(particles[dip->iAcol].dips.size()) != 1)
    return false;

  return true;
}

//--------------------------------------------------------------------------

// Check that trials do not contain junctions/ unusable pseudoparticles.

bool ColourReconnection::checkJunctionTrials() {
  for (int i = 0;i < int(junTrials.size());++i) {
    int minus = 0;
    if (junTrials[i].mode == 3)
      minus = 1;
    for (int j = 0;j < int(junTrials[i].dips.size()) - minus; ++j) {
      ColourDipole* dip = junTrials[i].dips[j];
      if (dip->isJun || dip->isAntiJun) {
        junTrials[i].list();
        return false;
      }
      if (particles[dip->iCol].dips.size() != 1 ||
          particles[dip->iAcol].dips.size() != 1) {
        junTrials[i].list();
        return false;
      }
    }
  }
  return true;
}


//--------------------------------------------------------------------------

// Find the neighbour to the colour side. Return false if the dipole
// is connected to a junction or the new particle has a junction inside of it.

bool ColourReconnection::findColNeighbour(ColourDipole*& dip) {
  // If only one active dipole, it has to be an antiquark.
  if (int(particles[dip->iCol].activeDips.size())  == 1)
    return false;

  // Has to have to active dipoles, otherwise something went wrong.
  if (int(particles[dip->iCol].activeDips.size())  != 2) {
    infoPtr->errorMsg("Warning in ColourReconnection::findAntiNeighbour: "
                      "Wrong number of active dipoles");
    return false;
  }
  // Otherwise find new dipole.
  if (dip == particles[dip->iCol].activeDips[0])
    dip = particles[dip->iCol].activeDips[1];
  else dip = particles[dip->iCol].activeDips[0];

  // Do not allow the new dipole to be connected to an antijunction.
  if (dip->isJun || dip->isAntiJun)
    return false;

  // Do not allow new dipole to have a pseudoparticle with
  // a baryon number inside.
  if (int(particles[dip->iCol].dips.size()) != 1)
    return false;

  return true;
}

//--------------------------------------------------------------------------

// Store used dipoles for a junction formation.

void ColourReconnection::storeUsedDips(TrialReconnection& trial) {
  // Normal dipole swap.
  if (trial.mode == 5) {

    for (int i = 0;i < 2; ++i) {
      ColourDipole* dip = trial.dips[i];
      if (dip->iCol < 0)
        for (int j = 0;j < 3; ++j)
        usedDipoles.push_back(junctions[-(dip->iCol / 10 + 1)].getColDip(j));
      if (dip->iAcol < 0)
        for (int j = 0;j < 3; ++j)
        usedDipoles.push_back(junctions[-(dip->iAcol / 10 + 1)].getColDip(j));

      usedDipoles.push_back(dip);
    }

  } else {

    for (int i = 0;i < 4; ++i) {
      if (trial.mode == 3 && i == 3)
        continue;
      usedDipoles.push_back(trial.dips[i]);
      ColourDipole* dip = trial.dips[i];


      while (findAntiNeighbour(dip) && dip != trial.dips[i])
        usedDipoles.push_back(dip);

      dip = trial.dips[i];
      while (findColNeighbour(dip) && dip != trial.dips[i])
        usedDipoles.push_back(dip);
    }
  }
}

//--------------------------------------------------------------------------

// Calculate the difference between the old and new lambda for dipole swap.

double ColourReconnection::getLambdaDiff(ColourDipole* dip1,
  ColourDipole* dip2) {

  // Needed to make sure the same dipoles are compared.
  vector<ColourDipole*> oldDips, newDips;

  // Calculate old string length.
  double oldLambda = calculateStringLength(dip1, oldDips)
    + calculateStringLength( dip2, oldDips);

  // Make test configuration.
  swapDipoles(dip1,dip2);

 // Calculate new string lengths
  double newLambda = calculateStringLength(dip1, newDips)
    +  calculateStringLength(dip2, newDips);

  // Swap back.
  swapDipoles(dip1, dip2, true);

  // First check if new combination was not useable.
  if (newLambda >= 0.5E9) return -1e9;

  // Return the difference.
  return oldLambda - newLambda;
}

//--------------------------------------------------------------------------

// Calculate the difference between the old and new lambda.

double ColourReconnection::getLambdaDiff(ColourDipole* dip1,
  ColourDipole* dip2, ColourDipole* dip3, ColourDipole* dip4, int mode) {

  // Calculate old lambda measure.

  double oldLambda = calculateStringLength(dip1->iCol, dip1->iAcol)
    + calculateStringLength(dip2->iCol, dip2->iAcol);
  if (dip3 != dip1) 
    oldLambda += calculateStringLength(dip3->iCol, dip3->iAcol);
  if (dip4 != 0 && dip4 != dip2)
    oldLambda += calculateStringLength(dip4->iCol, dip4->iAcol);

  // Calculate new lambda.
  double newLambda = 0;

  if (mode == 0)
      newLambda = calculateDoubleJunctionLength(dip1->iCol, dip2->iCol,
                                                dip1->iAcol, dip2->iAcol);
  else if (mode == 1) {
    if (dip2 == dip4)
      newLambda = calculateJunctionLength(dip1->iCol, dip2->iCol, dip3->iCol)
        + calculateJunctionLength(dip1->iAcol, dip2->iAcol, dip3->iAcol);
    else
      newLambda = calculateJunctionLength(dip1->iCol, dip2->iCol, dip3->iCol)
        + calculateJunctionLength(dip2->iAcol, dip3->iAcol, dip4->iAcol)
          + calculateStringLength(dip4->iCol, dip1->iAcol);
  }

  else if (mode == 2) {
    if (dip1 == dip3)
      newLambda = calculateJunctionLength(dip1->iCol, dip2->iCol, dip4->iCol)
        + calculateJunctionLength(dip1->iAcol, dip2->iAcol, dip4->iAcol);
    else
      newLambda = calculateJunctionLength(dip1->iCol, dip2->iCol, dip4->iCol)
          + calculateJunctionLength(dip1->iAcol, dip3->iAcol, dip4->iAcol)
          + calculateStringLength(dip3->iCol, dip2->iAcol);
  }

  // Triple junction connection.
  else if (mode == 3)
    newLambda = calculateJunctionLength(dip1->iCol, dip2->iCol, dip3->iCol)
      + calculateJunctionLength(dip1->iAcol, dip2->iAcol, dip3->iAcol);

  // First check if new combination was not useable.
  if (newLambda >= 0.5E9) return -1e9;

  // Returning result.
  return oldLambda - newLambda;

}

//--------------------------------------------------------------------------

// Change colour structure to describe the reconnection in juncTrial.

void ColourReconnection::doDipoleTrial(TrialReconnection& trial) {

  // Store for easier use.
  ColourDipole* dip1 = trial.dips[0];
  ColourDipole* dip2 = trial.dips[1];

  // If both acols ends are normal particles.
  if (dip1->iAcol >= 0 && dip2->iAcol >= 0) {
    swap(particles[dip1->iAcol].dips[dip1->iAcolLeg].front()->iAcol,
         particles[dip2->iAcol].dips[dip2->iAcolLeg].front()->iAcol);
    swap(particles[dip1->iAcol].dips[dip1->iAcolLeg].front(),
         particles[dip2->iAcol].dips[dip2->iAcolLeg].front());
  // If only dip1 has normal acol end.
  } else if (dip1->iAcol >= 0) {
    swap(particles[dip1->iAcol].dips[dip1->iAcolLeg].front()->iAcol,
       junctions[-(dip2->iAcol / 10 + 1)].dipsOrig[-dip2->iAcol % 10]->iAcol);
    swap(particles[dip1->iAcol].dips[dip1->iAcolLeg].front(),
         junctions[-(dip2->iAcol / 10 + 1)].dipsOrig[-dip2->iAcol % 10]);
  // If only dip2 has normal acol end.
  } else if (dip2->iAcol >= 0) {
    swap(particles[dip2->iAcol].dips[dip2->iAcolLeg].front()->iAcol,
       junctions[-(dip1->iAcol / 10 + 1)].dipsOrig[-dip1->iAcol % 10]->iAcol);
    swap(particles[dip2->iAcol].dips[dip2->iAcolLeg].front(),
         junctions[-(dip1->iAcol / 10 + 1)].dipsOrig[-dip1->iAcol % 10]);
  // If both ends are junctions.
  } else {
    swap(junctions[ -(dip1->iAcol / 10 + 1) ].dipsOrig[
           -dip1->iAcol % 10 ]->iAcol,
         junctions[ -(dip2->iAcol / 10 + 1) ].dipsOrig[
           -dip2->iAcol % 10 ]->iAcol);
    swap(junctions[ -(dip1->iAcol / 10 + 1) ].dipsOrig[ -dip1->iAcol % 10],
         junctions[ -(dip2->iAcol / 10 + 1) ].dipsOrig[ -dip2->iAcol % 10] );
  }

  // Swap the dipoles.
  swapDipoles(dip1, dip2);

  // If new particles are below treshhold, form pseudoParticles.
  if (mDip(dip1) < m0) makePseudoParticle(dip1, 110, true);
  if (mDip(dip2) < m0) makePseudoParticle(dip2, 110, true);

  // Done.

}

//--------------------------------------------------------------------------

// Update the list of dipole trial swaps to account for latest swap.

void ColourReconnection::updateDipoleTrials() {

  // Remove any dipTrials that contains an used dipole.
  for (int i = 0; i < int(dipTrials.size()); ++i)
    for (int j = 0;j < 2; ++j) {
      if (binary_search(usedDipoles.begin(), usedDipoles.end(),
                       dipTrials[i].dips[j]) ) {
        dipTrials.erase(dipTrials.begin() + i);
        i--;
        break;
      }
    }

  // Make list of active dipoles.
  vector<ColourDipole*> activeDipoles;
  for (int i = 0;i < int(dipoles.size()); ++i)
    if (dipoles[i]->isActive)
      activeDipoles.push_back(dipoles[i]);

  // Loop over list of used dipoles and create new trial reconnections.
  for (int i = 0;i < int(usedDipoles.size()); ++i)
    if (usedDipoles[i]->isActive)
      for (int j = 0; j < int(activeDipoles.size()); ++j)
        singleReconnection(usedDipoles[i], activeDipoles[j]);

}

//--------------------------------------------------------------------------

// Update the list of dipole trial swaps to account for latest swap.

void ColourReconnection::updateJunctionTrials() {

 // Remove any junTrials that contains an used dipole.
  for (int i = 0; i < int(junTrials.size()); ++i)
    for (int j = 0; j < 4; ++j) {
      if (binary_search(usedDipoles.begin(), usedDipoles.end(),
                       junTrials[i].dips[j]) ) {
        junTrials.erase(junTrials.begin() + i);
        i--;
        break;
      }
    }

  // Make list of active dipoles.
  vector<ColourDipole*> activeDipoles;
  for (int i = 0;i < int(dipoles.size()); ++i)
    if (dipoles[i]->isActive)
      activeDipoles.push_back(dipoles[i]);

  // Loop over used dipoles and form new junction trials.
  for (int i = 0;i < int(usedDipoles.size()); ++i)
    if (usedDipoles[i]->isActive)
      for (int j = 0; j < int(activeDipoles.size()); ++j)
        singleJunction(usedDipoles[i], activeDipoles[j]);

  // Loop over used dipoles and form new junction trials.
  for (int i = 0;i < int(usedDipoles.size()); ++i)
    if (usedDipoles[i]->isActive)
      for (int j = 0; j < int(activeDipoles.size()); ++j)
        for (int k = j + 1; k < int(activeDipoles.size()); ++k)
          singleJunction(usedDipoles[i], activeDipoles[j], activeDipoles[k]);

}

//--------------------------------------------------------------------------
        
// Change colour structure to describe the reconnection in juncTrial.
void ColourReconnection::doJunctionTrial(Event& event,
  TrialReconnection& juncTrial) {

  int mode = juncTrial.mode;
  // If trial mode is 3 (three dipoles -> 2 junctions) use its own update.
  if (mode == 3) {
    doTripleJunctionTrial(event, juncTrial);
    return;
  }

  // Store dipoles and numbers for easier acces.
  ColourDipole* dip1 = juncTrial.dips[0];
  ColourDipole* dip2 = juncTrial.dips[1];
  ColourDipole* dip3 = juncTrial.dips[2];
  ColourDipole* dip4 = juncTrial.dips[3];

  int iCol1 = dip1->iCol;
  int iCol2 = dip2->iCol;
  int iCol3 = dip3->iCol;
  int iCol4 = dip4->iCol;
  int iAcol1 = dip1->iAcol;
  int iAcol2 = dip2->iAcol;
  int iAcol3 = dip3->iAcol;
  int iAcol4 = dip4->iAcol;

  int oldCol1 = dip1->col;
  int oldCol2 = dip2->col;
  int oldCol3 = dip3->col;
  int oldCol4 = dip4->col;

  // New colour tags needed, since three more dipoles will be made.
  int newCol1 = event.nextColTag();
  int newCol2 = event.nextColTag();
  int newCol3 = event.nextColTag();

  // Need to make 3 new real dipoles and 3 active dipoles.

  // First make dipoles between junctions.

  // Find the junction colour.
  int junCol = 3 * (3 - (dip1->colReconnection / 3)
             - (dip2->colReconnection / 3) ) + dip1->colReconnection % 3;

  // if other than 9 colours.
  if (nReconCols != 9) {
    while (junCol < 0 || junCol % 3 != dip1->colReconnection % 3 ||
           junCol == dip1->colReconnection || junCol == dip2->colReconnection)
      junCol = int(nReconCols * rndmPtr->flat());
  }
  // Need one active and one real dipole.
  int iJun = junctions.size();
  int iAntiJun = junctions.size() + 1;

  // Store real dipoles.
  ColourDipole* dip3real = particles[iCol3].dips[dip3->iColLeg].back();
  ColourDipole* dip4real = particles[iCol4].dips[dip4->iColLeg].back();

  // If the junction and antijunction are directly connected.
  int iActive1 = 0, iReal1 = 0;
  if (mode == 0) {
    dipoles.push_back(new ColourDipole(newCol1, -( iAntiJun * 10 + 10 + 2) ,
      -( iJun * 10 + 10 + 2), junCol, true, true, false, true));
    iReal1 = dipoles.size() - 1;
    dipoles.push_back(new ColourDipole(newCol1, -( iAntiJun * 10 + 10 + 2) ,
      -( iJun * 10 + 10 + 2), junCol, true, true));
    iActive1 = dipoles.size() - 1;
  } else if (mode == 1) {
    int iCol3real = particles[iCol3].dips[dip3->iColLeg].back()->iCol;
     dipoles.push_back(new ColourDipole(newCol1, iCol3real ,
      -( iJun * 10 + 10 + 2), junCol, true, false, false, true));
    iReal1 = dipoles.size() - 1;
    particles[iCol3].dips[dip3->iColLeg].back() = dipoles.back();
    dipoles.push_back(new ColourDipole(newCol1, dip3->iCol,
      -( iJun * 10 + 10 + 2), junCol, true, false));
    iActive1 = dipoles.size() - 1;
  } else if (mode == 2) {
    int iCol4real = particles[iCol4].dips[dip4->iColLeg].back()->iCol;
    dipoles.push_back(new ColourDipole(newCol1, iCol4real,
      -( iJun * 10 + 10 + 2), junCol, true, false, false, true));
    iReal1 = dipoles.size() - 1;
    particles[iCol4].dips[dip4->iColLeg].back() = dipoles.back();
    dipoles.push_back(new ColourDipole(newCol1, dip4->iCol,
      -( iJun * 10 + 10 + 2), junCol, true, false));
    iActive1 = dipoles.size() - 1;
  }

  // Now make dipole between anti junction and iAcol1.
  // Start by finding real iAcol.
  int iAcol3real  = particles[iAcol3].dips[dip3->iAcolLeg].front()->iAcol;
  dipoles.push_back(new ColourDipole(newCol2, -( iAntiJun * 10 + 10),
    iAcol3real, dip3->colReconnection, false, true, false, true));
  int iReal2 = dipoles.size() - 1;
  particles[iAcol3].dips[dip3->iAcolLeg].front() = dipoles.back();

  dipoles.push_back(new ColourDipole(newCol2, -( iAntiJun * 10 + 10),
    iAcol3, dip3->colReconnection, false, true));
  dipoles.back()->iAcolLeg = dip3->iAcolLeg;
  int iActive2 = dipoles.size() - 1;

  // Now make dipole between anti junction and iAcol1.
  // Start by finding real iAcol.
  int iAcol4real = particles[iAcol4].dips[dip4->iAcolLeg].front()->iAcol;
  dipoles.push_back(new ColourDipole(newCol3, -( iAntiJun * 10 + 10 + 1),
    iAcol4real, dip4->colReconnection, false, true, false, true));
  int iReal3 = dipoles.size() - 1;
  particles[iAcol4].dips[dip4->iAcolLeg].front() = dipoles.back();

  dipoles.push_back(new ColourDipole(newCol3, -( iAntiJun * 10 + 10 + 1),
    iAcol4, dip4->colReconnection, false, true));
  dipoles.back()->iAcolLeg = dip4->iAcolLeg;
  int iActive3 = dipoles.size() - 1;

  // Update already existing dipoles, start by internal dipoles.
  // Now take dipoles connected to the anti junction
  // and a possible gluon-gluon connection.
  if (mode == 1) {
    if (dip2 == dip4) {

      // Update real dipole.
      dip3real->iAcol = particles[dip1->iAcol].dips[dip1->iAcolLeg].
        front()->iAcol;
      dip3real->iCol  = -( iAntiJun * 10 + 10 + 2);
      dip3real->isAntiJun = true;

      // Update active dipoles.
      dip3->iAcol = dip1->iAcol;
      dip3->iAcolLeg = dip1->iAcolLeg;
      dip3->isAntiJun = true;
      dip3->iCol = -( iAntiJun * 10 + 10 + 2);
      dip3->iColLeg = 0;

      // Store real dipole
      particles[dip3->iAcol].dips[dip3->iAcolLeg].front() = dip3real;

    } else {

      // Update real dipole.
      dip3real->iAcol = particles[dip2->iAcol].dips[dip2->iAcolLeg].
        front()->iAcol;
      dip3real->iCol  = -( iAntiJun * 10 + 10 + 2);
      dip3real->isAntiJun = true;
      dip4real->iAcol = particles[dip1->iAcol].dips[dip1->iAcolLeg].
        front()->iAcol;

      // Change the dipole connected to the antijunction.
      dip3->iAcol = dip2->iAcol;
      dip3->iAcolLeg = dip2->iAcolLeg;
      dip3->isAntiJun = true;
      dip3->iCol = -( iAntiJun * 10 + 10 + 2);
      dip3->iColLeg = 0;

      // Change the dipole between the two gluons.
      dip4->iAcol = dip1->iAcol;
      dip4->iAcolLeg = dip1->iAcolLeg;

      // Store real dipole
      particles[dip3->iAcol].dips[dip3->iAcolLeg].front() = dip3real;
      particles[dip4->iAcol].dips[dip4->iAcolLeg].front() = dip4real;

    }
  } else if (mode == 2) {
    if (dip1 == dip3) {

      // Update real dipole.
      dip4real->iAcol = particles[dip2->iAcol].dips[dip2->iAcolLeg].
        front()->iAcol;
      dip4real->iCol  = -( iAntiJun * 10 + 10 + 2);
      dip4real->isAntiJun = true;

      // Update active dipoles.
      dip4->iAcol = dip2->iAcol;
      dip4->iAcolLeg = dip2->iAcolLeg;
      dip4->isAntiJun = true;
      dip4->iCol = -( iAntiJun * 10 + 10 + 2);
      dip4->iColLeg = 0;

      // Store real dipole
      particles[dip4->iAcol].dips[dip4->iAcolLeg].front() = dip4real;

    } else {

      // Update real dipole.
      dip4real->iAcol = particles[dip1->iAcol].dips[dip1->iAcolLeg].
        front()->iAcol;
      dip4real->iCol  = -( iAntiJun * 10 + 10 + 2);
      dip4real->isAntiJun = true;
      dip3real->iAcol = particles[dip2->iAcol].dips[dip2->iAcolLeg].
        front()->iAcol;

      // Change the dipole connected to the antijunction.
      dip4->iAcol = dip1->iAcol;
      dip4->iAcolLeg = dip1->iAcolLeg;
      dip4->isAntiJun = true;
      dip4->iCol = -( iAntiJun * 10 + 10 + 2);
      dip4->iColLeg = 0;

      // Change the dipole between the two gluons.
      dip3->iAcol = dip2->iAcol;
      dip3->iAcolLeg = dip2->iAcolLeg;

      // Store real dipole
      particles[dip3->iAcol].dips[dip3->iAcolLeg].front() = dip3real;
      particles[dip4->iAcol].dips[dip4->iAcolLeg].front() = dip4real;
    }
  }

  // Dipoles connected to the junction.
  // Update real dipoles.
  particles[iCol1].dips[dip1->iColLeg].back()->iAcol = - (iJun * 10 + 10);
  particles[iCol2].dips[dip2->iColLeg].back()->iAcol = - (iJun * 10 + 10 + 1);
  particles[iCol1].dips[dip1->iColLeg].back()->isJun = true;
  particles[iCol2].dips[dip2->iColLeg].back()->isJun = true;

  // Update active dipoles.
  dip1->isJun = true;
  dip2->isJun = true;
  dip1->iAcol = - (iJun * 10 + 10);
  dip2->iAcol = - (iJun * 10 + 10 + 1);
  dip1->iAcolLeg = 0;
  dip2->iAcolLeg = 0;

  // Update active dipoles for anti particles.
  // Normally should only contain active dipoles once,
  // only problem is if the two dipole ends are the same particle.
  // Start by settings common dipoles.
  for (int i = 0; i < int(particles[iAcol3].activeDips.size()); ++i)
    if (particles[iAcol3].activeDips[i] == dip3) {
      particles[iAcol3].activeDips[i] = dipoles[iActive2];
      break;
    }

  for (int i = 0; i < int(particles[iAcol4].activeDips.size()); ++i)
    if (particles[iAcol4].activeDips[i] == dip4) {
      particles[iAcol4].activeDips[i] = dipoles[iActive3];
    break;
    }

  // Depending on how the new string is connected, the active dipoles vary.
  if (mode == 1) {
    for (int i = 0; i < int(particles[iCol3].activeDips.size()); ++i)
      if (particles[iCol3].activeDips[i] == dip3) {
        particles[iCol3].activeDips[i] = dipoles[iActive1];
        break;
      }

    if (dip2 == dip4) {
      for (int i = 0; i < int(particles[iAcol1].activeDips.size()); ++i)
        if (particles[iAcol1].activeDips[i] == dip1) {
          particles[iAcol1].activeDips[i] = dip3;
          break;
        }
    } else {
      for (int i = 0; i < int(particles[iAcol2].activeDips.size()); ++i)
        if (particles[iAcol2].activeDips[i] == dip2) {
          particles[iAcol2].activeDips[i] = dip3;
          break;
        }

      for (int i = 0; i < int(particles[iAcol1].activeDips.size()); ++i)
        if (particles[iAcol1].activeDips[i] == dip1) {
          particles[iAcol1].activeDips[i] = dip4;
          break;
        }
    }
  } else if (mode == 2) {
    for (int i = 0; i < int(particles[iCol4].activeDips.size()); ++i)
      if (particles[iCol4].activeDips[i] == dip4) {
        particles[iCol4].activeDips[i] = dipoles[iActive1];
        break;
      }

    if (dip1 == dip3) {
      for (int i = 0; i < int(particles[iAcol2].activeDips.size()); ++i)
        if (particles[iAcol2].activeDips[i] == dip2) {
          particles[iAcol2].activeDips[i] = dip4;
          break;
        }
    } else {
      for (int i = 0; i < int(particles[iAcol1].activeDips.size()); ++i)
        if (particles[iAcol1].activeDips[i] == dip1) {
          particles[iAcol1].activeDips[i] = dip4;
          break;
        }

      for (int i = 0; i < int(particles[iAcol2].activeDips.size()); ++i)
        if (particles[iAcol2].activeDips[i] == dip2) {
          particles[iAcol2].activeDips[i] = dip3;
          break;
        }
    }
  }

  // Add the junctions to the event.
  junctions.push_back(Junction(1, oldCol1, oldCol2, newCol1));
  if (mode == 0) junctions.push_back(Junction(2, newCol2, newCol3, newCol1));
  else if (mode == 1)
    junctions.push_back(Junction(2, newCol2, newCol3, oldCol3));
  else if (mode == 2)
    junctions.push_back(Junction(2, newCol2, newCol3, oldCol4));

  // Set junction information.
  junctions[iJun].dipsOrig[0] =
    particles[iCol1].dips[dip1->iColLeg].back();
  junctions[iJun].dipsOrig[1] =
    particles[iCol2].dips[dip2->iColLeg].back();
  junctions[iJun].dipsOrig[2] = dipoles[iReal1];
  junctions[iJun].dips[0] = dip1;
  junctions[iJun].dips[1] = dip2;
  junctions[iJun].dips[2] = dipoles[iActive1];

  // Set anti junction information.
  junctions[iAntiJun].dips[0] = dipoles[iActive2];
  junctions[iAntiJun].dips[1] = dipoles[iActive3];
  junctions[iAntiJun].dipsOrig[0] = dipoles[iReal2];
  junctions[iAntiJun].dipsOrig[1] = dipoles[iReal3];

  if (mode == 0) {
    junctions[iAntiJun].dips[2] = dipoles[iActive1];
    junctions[iAntiJun].dipsOrig[2] = dipoles[iReal1];
  } else if (mode == 1) {
    junctions[iAntiJun].dips[2] = dip3;
    junctions[iAntiJun].dipsOrig[2] =
      particles[dip3->iAcol].dips[dip3->iAcolLeg].front();
  } else if (mode == 2) {
    junctions[iAntiJun].dips[2] = dip4;
    junctions[iAntiJun].dipsOrig[2] =
      particles[dip4->iAcol].dips[dip4->iAcolLeg].front();
  }

  // Make pseudo particles.
  if (dip1->isActive && mDip(dip1) < m0)
    makePseudoParticle(dip1, 110, true);
  if (dip2->isActive && mDip(dip2) < m0)
    makePseudoParticle(dip2, 110, true);
  if (dip3->isActive && mDip(dip3) < m0)
    makePseudoParticle(dip3, 110, true);
  if (dip4->isActive && mDip(dip4) < m0)
    makePseudoParticle(dip4, 110, true);

  if (dipoles[iActive1]->isActive && mDip(dipoles[iActive1]) < m0)
    makePseudoParticle(dipoles[iActive1], 110, true);
  if (dipoles[iActive2]->isActive && mDip(dipoles[iActive2]) < m0)
    makePseudoParticle(dipoles[iActive2], 110, true);
  if (dipoles[iActive3]->isActive && mDip(dipoles[iActive3]) < m0)
    makePseudoParticle(dipoles[iActive3], 110, true);

  // Add new dipoles to usedDipoles.
  usedDipoles.push_back(dipoles[iActive1]);
  usedDipoles.push_back(dipoles[iActive2]);
  usedDipoles.push_back(dipoles[iActive3]);

  // Done.
}

//--------------------------------------------------------------------------

void ColourReconnection::doTripleJunctionTrial(Event& event,
  TrialReconnection& juncTrial) {

  // store information for easier acces.
  ColourDipole* dip1 = juncTrial.dips[0];
  ColourDipole* dip2 = juncTrial.dips[1];
  ColourDipole* dip3 = juncTrial.dips[2];

  // Store indices.
  int iCol1 = dip1->iCol;
  int iCol2 = dip2->iCol;
  int iCol3 = dip3->iCol;
  int iAcol1 = dip1->iAcol;
  int iAcol2 = dip2->iAcol;
  int iAcol3 = dip3->iAcol;

  // Store colours
  int oldCol1 = dip1->col;
  int oldCol2 = dip2->col;
  int oldCol3 = dip3->col;
  int newCol1 = event.nextColTag();
  int newCol2 = event.nextColTag();
  int newCol3 = event.nextColTag();

  // Store new junction indices.
  int iJun = junctions.size();
  int iAntiJun = junctions.size() + 1;

  // Now make dipole between anti junction and iAcol1.
  // Start by finding real iAcol.
  int iAcol1real
    = particles[iAcol1].dips[dip1->iAcolLeg].front()->iAcol;
  dipoles.push_back(new ColourDipole(newCol1, -( iAntiJun * 10 + 10),
    iAcol1real, dip1->colReconnection, false, true, false, true));
  int iReal1 = dipoles.size() - 1;
  particles[iAcol1].dips[dip1->iAcolLeg].front() = dipoles.back();

  dipoles.push_back(new ColourDipole(newCol1, -( iAntiJun * 10 + 10),
    iAcol1, dip1->colReconnection, false, true));
  dipoles.back()->iAcolLeg = dip1->iAcolLeg;
  int iActive1 = dipoles.size() - 1;

  // Now make dipole between anti junction and iAcol2.
  // Start by finding real iAcol2.
  int iAcol2real
    = particles[iAcol2].dips[dip2->iAcolLeg].front()->iAcol;
  dipoles.push_back(new ColourDipole(newCol2, -( iAntiJun * 10 + 10 + 1),
    iAcol2real, dip2->colReconnection, false, true, false, true));
  int iReal2 = dipoles.size() - 1;
  particles[iAcol2].dips[dip2->iAcolLeg].front() = dipoles.back();

  dipoles.push_back(new ColourDipole(newCol2, -( iAntiJun * 10 + 10 + 1),
    iAcol2, dip2->colReconnection, false, true));
  dipoles.back()->iAcolLeg = dip2->iAcolLeg;
  int iActive2 = dipoles.size() - 1;

  // Now make dipole between anti junction and iAcol3.
  // Start by finding real iAcol3.
  int iAcol3real
    = particles[iAcol3].dips[dip3->iAcolLeg].front()->iAcol;
  dipoles.push_back(new ColourDipole(newCol3, -( iAntiJun * 10 + 10 + 2),
    iAcol3real, dip3->colReconnection, false, true, false, true));
  int iReal3 = dipoles.size() - 1;
  particles[iAcol3].dips[dip3->iAcolLeg].front() = dipoles.back();

  dipoles.push_back(new ColourDipole(newCol3, -( iAntiJun * 10 + 10 + 2),
    iAcol3, dip3->colReconnection, false, true));
  dipoles.back()->iAcolLeg = dip3->iAcolLeg;
  int iActive3 = dipoles.size() - 1;

  // Update already existing dipoles.

  // Update real dipoles.
  particles[iCol1].dips[dip1->iColLeg].back()->iAcol = - (iJun * 10 + 10);
  particles[iCol2].dips[dip2->iColLeg].back()->iAcol = - (iJun * 10 + 10 + 1);
  particles[iCol3].dips[dip3->iColLeg].back()->iAcol = - (iJun * 10 + 10 + 2);
  particles[iCol1].dips[dip1->iColLeg].back()->isJun = true;
  particles[iCol2].dips[dip2->iColLeg].back()->isJun = true;
  particles[iCol3].dips[dip3->iColLeg].back()->isJun = true;

  // Update active dipoles.
  dip1->isJun = true;
  dip2->isJun = true;
  dip3->isJun = true;
  dip1->iAcol = - (iJun * 10 + 10);
  dip2->iAcol = - (iJun * 10 + 10 + 1);
  dip3->iAcol = - (iJun * 10 + 10 + 2);
  dip1->iAcolLeg = 0;
  dip2->iAcolLeg = 0;
  dip3->iAcolLeg = 0;

  // Update active dipoles for anti particles.
  for (int i = 0; i < int(particles[iAcol1].activeDips.size()); ++i)
    if (particles[iAcol1].activeDips[i] == dip1)
      particles[iAcol1].activeDips[i] = dipoles[iActive1];
  for (int i = 0; i < int(particles[iAcol2].activeDips.size()); ++i)
    if (particles[iAcol2].activeDips[i] == dip2)
      particles[iAcol2].activeDips[i] = dipoles[iActive2];
  for (int i = 0; i < int(particles[iAcol3].activeDips.size()); ++i)
    if (particles[iAcol3].activeDips[i] == dip3)
      particles[iAcol3].activeDips[i] = dipoles[iActive3];

  // Add the junctions to the event.
  junctions.push_back(Junction(1, oldCol1, oldCol2, oldCol3));
  junctions.push_back(Junction(2, newCol1, newCol3, newCol3));

  // Update junction ends.
  junctions[iJun].dipsOrig[0] =
    particles[iCol1].dips[dip1->iColLeg].back();
  junctions[iJun].dipsOrig[1] =
    particles[iCol2].dips[dip2->iColLeg].back();
  junctions[iJun].dipsOrig[2] =
    particles[iCol3].dips[dip3->iColLeg].back();
  junctions[iJun].dips[0] = dip1;
  junctions[iJun].dips[1] = dip2;
  junctions[iJun].dips[2] = dip3;

  // Update the anti junction.
  junctions[iAntiJun].dips[0] = dipoles[iActive1];
  junctions[iAntiJun].dips[1] = dipoles[iActive2];
  junctions[iAntiJun].dips[2] = dipoles[iActive3];
  junctions[iAntiJun].dipsOrig[0] = dipoles[iReal1];
  junctions[iAntiJun].dipsOrig[1] = dipoles[iReal2];
  junctions[iAntiJun].dipsOrig[2] = dipoles[iReal3];

  // Make pseudo particles if needed.
  if (dip1->isActive && mDip(dip1) < m0)
    makePseudoParticle(dip1, 110, true);
  if (dip2->isActive && mDip(dip2) < m0)
    makePseudoParticle(dip2, 110, true);
  if (dip3->isActive && mDip(dip3) < m0)
    makePseudoParticle(dip3, 110, true);

  if (dipoles[iActive1]->isActive && mDip(dipoles[iActive1]) < m0)
    makePseudoParticle(dipoles[iActive1], 110, true);
  if (dipoles[iActive2]->isActive && mDip(dipoles[iActive2]) < m0)
    makePseudoParticle(dipoles[iActive2], 110, true);
  if (dipoles[iActive3]->isActive && mDip(dipoles[iActive3]) < m0)
    makePseudoParticle(dipoles[iActive3], 110, true);

  // Add to newly created dipoles to used dipoles.
  usedDipoles.push_back(dipoles[iActive1]);
  usedDipoles.push_back(dipoles[iActive2]);
  usedDipoles.push_back(dipoles[iActive3]);

  // Done.
}

//--------------------------------------------------------------------------

// Allow colour reconnections by moving gluons from their current location
// to another colour line. Also optionally flip two colour chains.

bool ColourReconnection::reconnectMove( Event&  event, int oldSize) {

  // Create or reset arrays to prepare for the new event analysis.
  vector<int> iGlu;
  iReduceCol.resize( event.size() );
  iExpandCol.clear();
  map<int, int> colMap, acolMap;
  map<int, int>::iterator colM, acolM;
  vector<InfoGluonMove> infoGM;

  // Temporary variables.
  int iNow            = 0;
  int colNow          = 0;
  int acolNow         = 0;
  int iColNow         = 0;
  int iAcolNow        = 0;
  int col2Now         = 0;
  int iCol2Now        = 0;
  int iAcol2Now       = 0;
  double lambdaRefNow = 0.;
  double dLambdaNow   = 0.;

  // Loop over all final particles. Store (fraction of) gluons to move.
  for (int i = oldSize; i < event.size(); ++i) if (event[i].isFinal()) {
    if (event[i].id() == 21 && rndmPtr->flat() < fracGluon)
      iGlu.push_back(i);

    // Store location of all colour and anticolour particles and indices.
    if (event[i].col() > 0 || event[i].acol() > 0) {
      iReduceCol[i] = iExpandCol.size();
      iExpandCol.push_back(i);
      if (event[i].col() > 0) colMap[event[i].col()] = i;
      if (event[i].acol() > 0) acolMap[event[i].acol()] = i;
    }
  }

  // Erase (anti)colours for (anti)junctions and skip adjacent gluons.
  for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {
    if (event.kindJunction(iJun) == 1) {
      for (int j = 0; j < 3; ++j) {
        int jCol = event.colJunction( iJun, j);
        for (colM = colMap.begin(); colM != colMap.end(); ++colM)
        if (colM->first == jCol) {
          colMap.erase( colM);
          break;
        }
        for (int iG = 0; iG < int(iGlu.size()); ++iG)
        if (event[iGlu[iG]].col() == jCol) {
          iGlu.erase(iGlu.begin() + iG);
          break;
        }
      }
    } else if (event.kindJunction(iJun) == 2) {
      for (int j = 0; j < 3; ++j) {
        int jCol = event.colJunction( iJun, j);
        for (acolM = acolMap.begin(); acolM != acolMap.end(); ++acolM)
        if (acolM->first == jCol) {
          acolMap.erase( acolM);
          break;
        }
        for (int iG = 0; iG < int(iGlu.size()); ++iG)
        if (event[iGlu[iG]].acol() == jCol) {
          iGlu.erase(iGlu.begin() + iG);
          break;
        }
      }
    }
  }

  // Error checks.
  int nGlu = iGlu.size();
  int nCol = colMap.size();
  if (int(acolMap.size()) != nCol) {
    infoPtr->errorMsg("Error in MBReconUserHooks: map sizes do not match");
    return false;
  }
  colM  = colMap.begin();
  acolM = acolMap.begin();
  for (int iCol = 0; iCol < nCol; ++iCol) {
    if (colM->first != acolM->first) {
      infoPtr->errorMsg("Error in MBReconUserHooks: map elements"
        " do not match");
      return false;
    }
    ++colM;
    ++acolM;
  }

  // Calculate and tabulate lambda between any pair of coloured partons.
  nColMove = iExpandCol.size();
  lambdaijMove.resize( pow2(nColMove) );
  for (int iAC = 0; iAC < nColMove - 1; ++iAC) {
    int i = iExpandCol[iAC];
    for (int jAC = iAC + 1; jAC < nColMove; ++jAC) {
      int j = iExpandCol[jAC];
      lambdaijMove[nColMove * iAC + jAC]
        = log(1. + m2( event[i], event[j]) / m2Lambda);
    }
  }

  // Set up initial possible gluon moves with lambda gains/losses.
  for (int iG = 0; iG < nGlu; ++iG) {

    // Gluon and its current neighbours.
    iNow     = iGlu[iG];
    colNow   = event[iNow].col();
    acolNow  = event[iNow].acol();
    iColNow  = acolMap[colNow];
    iAcolNow = colMap[acolNow];

    // Addition to Lambda of gluon in current position.
    lambdaRefNow = lambda123Move( iNow, iColNow, iAcolNow);

    // Loop over all colour lines where gluon could be inserted.
    for (colM = colMap.begin(); colM != colMap.end(); ++colM) {
      col2Now   = colM->first;
      iCol2Now  = colMap[col2Now];
      iAcol2Now = acolMap[col2Now];

      // Addition to total Lambda if gluon moved to be inserted on line.
      dLambdaNow = (iCol2Now == iNow || iAcol2Now == iNow
        || iColNow == iAcolNow) ? 2e4
        : lambda123Move( iNow, iCol2Now, iAcol2Now) - lambdaRefNow;

      // Add new container for gluon and colour line information.
      infoGM.push_back( InfoGluonMove( iNow, colNow, acolNow, iColNow,
        iAcolNow, col2Now, iCol2Now, iAcol2Now, lambdaRefNow, dLambdaNow ));
    }
  }
  int nPair = infoGM.size();

  // Keep on looping over moves until no further negative dLambda.
  for ( int iMove = 0; iMove < nGlu; ++iMove) {
    int    iPairMin   = -1;
    double dLambdaMin = 1e4;

    // Find lowest dLambda.
    for (int iPair = 0; iPair < nPair; ++iPair)
    if (infoGM[iPair].dLambda < dLambdaMin) {
      iPairMin   = iPair;
      dLambdaMin = infoGM[iPair].dLambda;
    }

    // Break if no shift below upper limit found.
    if (dLambdaMin > -dLambdaCut) break;

    // Partons and colours involved in move.
    InfoGluonMove& selSM = infoGM[iPairMin];
    int i1Sel     = selSM.i1;
    int iCol1Sel  = selSM.iCol1;
    int iAcol1Sel = selSM.iAcol1;
    int iCol2Sel  = selSM.iCol2;
    int iAcol2Sel = selSM.iAcol2;
    int iCol2Mod[3]  = { iAcol1Sel , i1Sel     , iCol2Sel    };
    int col2Mod[3]   = { selSM.col1, selSM.col2, selSM.acol1};

    // Remove gluon from old colour line and insert on new colour line.
    for (int i = 0; i < 3; ++i) {
      event[ iCol2Mod[i] ].col( col2Mod[i] );
      colMap[ col2Mod[i] ] = iCol2Mod[i];
    }

    // Update info for partons with new colors.
    int  i1Now    = 0;
    bool doUpdate = false;
    for (int iPair = 0; iPair < nPair; ++iPair) {
      InfoGluonMove& tmpSM = infoGM[iPair];
      if (tmpSM.i1 != i1Now) {
        i1Now = tmpSM.i1;
        doUpdate = false;
        if (i1Now == i1Sel || i1Now == iCol1Sel || i1Now == iAcol1Sel
          || i1Now == iCol2Sel || i1Now == iAcol2Sel) {
          colNow       = event[i1Now].col();
          acolNow      = event[i1Now].acol();
          iColNow      = acolMap[colNow];
          iAcolNow     = colMap[acolNow];
          lambdaRefNow = lambda123Move( i1Now, iColNow, iAcolNow);
          doUpdate     = true;
        }
      }
      if (doUpdate) {
        tmpSM.col1      = colNow;
        tmpSM.acol1     = acolNow;
        tmpSM.iCol1     = iColNow;
        tmpSM.iAcol1    = iAcolNow;
        tmpSM.lambdaRef = lambdaRefNow;
      }
    }

    // Update info on dLambda for affected particles and colour lines.
    for (int iPair = 0; iPair < nPair; ++iPair) {
      InfoGluonMove& tmpSM = infoGM[iPair];
      int iMod = -1;
      for (int i = 0; i < 3; ++i) if (tmpSM.col2 == col2Mod[i]) iMod = i;
      if (iMod > -1) tmpSM.iCol2 = iCol2Mod[iMod];
      if (tmpSM.i1 == i1Sel || tmpSM.i1 == iCol1Sel || tmpSM.i1 == iAcol1Sel
        || tmpSM.i1 == iCol2Sel || tmpSM.i1 == iAcol2Sel || iMod > -1)
        tmpSM.dLambda = (tmpSM.iCol2 == tmpSM.i1 || tmpSM.iAcol2 == tmpSM.i1
          || tmpSM.iCol1 == tmpSM.iAcol1) ? 2e4
          : lambda123Move( tmpSM.i1, tmpSM.iCol2, tmpSM.iAcol2)
          - tmpSM.lambdaRef;
    }

  // End of loop over gluon shifting.
  }

  // Done if no flip.
  if (flipMode == 0) return true;

  // Array with colour lines, and where each line begins and ends.
  vector<int> iTmpFlip, iVecFlip, iBegFlip, iEndFlip;

  // Variables for minimum search.
  int i1c, i1a, i2c, i2a, i1cMin, i1aMin, i2cMin, i2aMin, iSMin;
  double dLambdaFlip, dLambdaFlipMin;
  vector<InfoGluonMove> flipMin;

  // Grab all colour ends.
  for (int i = oldSize; i < event.size(); ++i)
  if (event[i].isFinal() && event[i].col() > 0 && event[i].acol() == 0) {
    iTmpFlip.clear();
    iTmpFlip.push_back( i);

    // Step through colour neighbours to catch system.
    iNow = i;
    acolM = acolMap.find( event[iNow].col() );
    bool foundEnd = false;
    while (acolM != acolMap.end()) {
      iNow = acolM->second;
      iTmpFlip.push_back( iNow);
      if (event[iNow].col() == 0) {
        foundEnd = true;
        break;
      }
      acolM = acolMap.find( event[iNow].col() );
    }

    // Store acceptable system, optionally including junction legs.
    if (foundEnd || flipMode == 2) {
      iBegFlip.push_back( iVecFlip.size());
      for (int j = 0; j < int(iTmpFlip.size()); ++j)
        iVecFlip.push_back( iTmpFlip[j]);
      iEndFlip.push_back( iVecFlip.size());
    }
  }

  // Optionally search for antijunction legs: grab all anticolour ends.
  if (flipMode == 2) for (int i = oldSize; i < event.size(); ++i)
  if (event[i].isFinal() && event[i].acol() > 0 && event[i].col() == 0) {
    iTmpFlip.clear();
    iTmpFlip.push_back( i);

    // Step through anticolour neighbours to catch system.
    iNow = i;
    colM = colMap.find( event[iNow].acol() );
    bool foundEnd = false;
    while (colM != colMap.end()) {
      iNow = colM->second;
      iTmpFlip.push_back( iNow);
      if (event[iNow].acol() == 0) {
        foundEnd = true;
        break;
      }
      colM = colMap.find( event[iNow].acol() );
    }

    // Store acceptable system, but do not doublecount q - (n g) - qbar.
    if (!foundEnd) {
      iBegFlip.push_back( iVecFlip.size());
      for (int j = 0; j < int(iTmpFlip.size()); ++j)
        iVecFlip.push_back( iTmpFlip[j]);
      iEndFlip.push_back( iVecFlip.size());
    }
  }

  // Loop through all system pairs.
  int nSysFlip = iBegFlip.size();
  for (int iSys1 = 0; iSys1 < nSysFlip - 1; ++iSys1)
  if (iBegFlip[iSys1] >= 0)
  for (int iSys2 = iSys1 + 1; iSys2 < nSysFlip; ++iSys2)
  if (iBegFlip[iSys2] >= 0) {
    i1cMin     = 0;
    i1aMin     = 0;
    i2cMin     = 0;
    i2aMin     = 0;
    dLambdaFlipMin = 1e4;

    // Loop through all possible flip locations for a pair.
    for (int j1 = iBegFlip[iSys1]; j1 < iEndFlip[iSys1] - 1; ++j1)
    for (int j2 = iBegFlip[iSys2]; j2 < iEndFlip[iSys2] - 1; ++j2) {
      i1c = iVecFlip[j1];
      i1a = iVecFlip[j1 + 1];
      i2c = iVecFlip[j2];
      i2a = iVecFlip[j2 + 1];
      dLambdaFlip = lambda12Move( i1c, i2a) + lambda12Move( i2c, i1a)
                  - lambda12Move( i1c, i1a) - lambda12Move( i2c, i2a);
      if (dLambdaFlip < dLambdaFlipMin) {
        i1cMin = i1c;
        i1aMin = i1a;
        i2cMin = i2c;
        i2aMin = i2a;
        dLambdaFlipMin = dLambdaFlip;
      }
    }

    // Store possible flips if low enough dLambdaMin.
    if (dLambdaFlipMin < -dLambdaCut) flipMin.push_back( InfoGluonMove(
      iSys1, iSys2, i1cMin, i1aMin, i2cMin, i2aMin, dLambdaFlipMin) );
  }
  int flipSize = flipMin.size();

  // Search for lowest possible flip among unused systems.
  for (int iFlip = 0; iFlip < min( nSysFlip / 2, flipSize); ++iFlip) {
    iSMin = -1;
    dLambdaFlipMin  = 1e4;
    for (int iSys12 = 0; iSys12 < flipSize; ++iSys12)
    if (flipMin[iSys12].i1 >= 0 && flipMin[iSys12].dLambda < dLambdaFlipMin) {
      iSMin   = iSys12;
      dLambdaFlipMin = flipMin[iSys12].dLambda;
    }

    // Do flip. Mark flipped systems.
    if (iSMin >= 0) {
      InfoGluonMove& flipNow = flipMin[iSMin];
      int iS1 = flipNow.i1;
      int iS2 = flipNow.i2;
      event[ flipNow.iAcol1 ].acol( event[flipNow.iCol2].col() );
      event[ flipNow.iAcol2 ].acol( event[flipNow.iCol1].col() );
       for (int iSys12 = 0; iSys12 < flipSize; ++iSys12)
      if ( flipMin[iSys12].i1 == iS1 || flipMin[iSys12].i1 == iS2
        || flipMin[iSys12].i2 == iS1 || flipMin[iSys12].i2 == iS2)
        flipMin[iSys12].i1 = -1;
    }
    else break;
  }

  // Done.
  return true;

}

//==========================================================================

} // end namespace Pythia8
