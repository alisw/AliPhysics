// ColourReconnectionHooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains two UserHooks that, along with the internal models,
// implement all the models used for the top mass study in
// S. Argyropoulos and T. Sjostrand,
// arXiv:1407.6653 [hep-ph] (LU TP 14-23, DESY 14-134, MCnet-14-15)

// MBReconUserHooks: can be used for all kinds of events, not only top ones.
// TopReconUserHooks: models intended specifically for top decay products,
// whereas the underlying event is handled by the default model.

// Warning: some small modifications have been made when collecting
// the models, but nothing intended to change the behaviour.
// Note: the move model is also available with ColourReconnection:mode = 2,
// while the ColourReconnection:mode = 1 model has not been used here.
// Note: the new models tend to be slower than the default CR scenario,
// since they have to probe many more reconnection possibilities.

#ifndef Pythia8_ColourReconnectionHooks_H
#define Pythia8_ColourReconnectionHooks_H

// Includes
#include "Pythia8/Pythia.h"
namespace Pythia8 {

//==========================================================================

// Class for colour reconnection models of general validity.

class MBReconUserHooks : public UserHooks {

public:

  // Constructor and destructor.
  // mode = 0: no reconnection (dummy option, does nothing);
  //      = 1: swap gluons to minimize lambda.
  //      = 2: move gluons to minimize lambda.
  // flip = 0: no flip between quark-antiquark ends.
  //      = 1: flip between quark-antiquark ends, excluding junction systems.
  //      = 2: flip between quark-antiquark ends, including junction systems.
  // dLamCut: smallest -delta-lambda value for which to swap/mode (positive).
  // fracGluon: the fraction of gluons that will be studied for reconnection.
  // m2Ref   : squared reference mass scale for lambda measure calculation.
  MBReconUserHooks(int modeIn = 0, int flipIn = 0, double dLamCutIn = 0.,
    double fracGluonIn = 1.) : mode(modeIn), flip(flipIn), dLamCut(dLamCutIn),
    fracGluon(fracGluonIn) { m2Ref = 1.; dLamCut = max(0., dLamCut); }
  ~MBReconUserHooks() {}

  // Allow colour reconnection after resonance decays (early or late)...
  virtual bool canReconnectResonanceSystems() {return true;}

  // ...which gives access to the event, for modification.
  virtual bool doReconnectResonanceSystems( int, Event& event) {

    // Return without action for relevant mode numbers.
    if (mode <= 0 || mode > 2) return true;

    // Double diffraction not yet implemented, so return without action.
    // (But works for internal move implementation.)
    if (infoPtr->isDiffractiveA() && infoPtr->isDiffractiveB()) return true;

    // Initial setup: relevant gluons and coloured partons.
    if (!setupConfig( event)) return false;

    // Done if not enough gluons.
    if ( (mode == 1 && nGlu < 2) || (mode == 2 && nGlu < 1) ) return true;

    // Colour reconnect. Return if failed.
    bool hasRec = (mode == 1) ? doReconnectSwap( event)
                              : doReconnectMove( event);
    if (!hasRec) return false;

    // Colour flip afterburner.
    if (flip > 0) return doReconnectFlip( event);
    return true;

  }

  // Return number of reconnections for current event.
  //int numberReconnections() {return nRec;}
  //double dLambdaReconnections() {return -dLamTot;}

private:

  // Mode. Number of reconnections. lambda measure reference scale.
  int    mode, flip, nRec, nGlu, nAllCol, nCol;
  double dLamCut, fracGluon, m2Ref, dLamTot;

  // Array of (indices of) final gluons.
  vector<int> iGlu;

  // Array of (indices of) all final coloured particles.
  vector<int> iToAllCol, iAllCol;

  // Maps telling where all colours and anticolours are stored.
  map<int, int> colMap, acolMap;

  // Array of all lambda distances between coloured partons.
  vector<double> lambdaij;

  // Function to return lambda value from array.
  double lambda12( int i, int j) {
    int iAC = iToAllCol[i]; int jAC = iToAllCol[j];
    return lambdaij[nAllCol * min( iAC, jAC) + max( iAC, jAC)];
  }

  // Function to return lambda(i,j) + lambda(i,k) - lambda(j,k).
  double lambda123( int i, int j, int k) {
    int iAC = iToAllCol[i]; int jAC = iToAllCol[j]; int kAC = iToAllCol[k];
    return lambdaij[nAllCol * min( iAC, jAC) + max( iAC, jAC)]
         + lambdaij[nAllCol * min( iAC, kAC) + max( iAC, kAC)]
         - lambdaij[nAllCol * min( jAC, kAC) + max( jAC, kAC)];
  }

  // Small nested class for the effect of a potential gluon swap or move.
  class InfoSwapMove{
    public:
    InfoSwapMove(int i1in = 0, int i2in = 0) : i1(i1in), i2(i2in) {}
    InfoSwapMove(int i1in, int i2in, int iCol1in, int iAcol1in, int iCol2in,
      int iAcol2in, int dLamIn) : i1(i1in), i2(i2in), iCol1(iCol1in),
      iAcol1(iAcol1in), iCol2(iCol2in), iAcol2(iAcol2in), dLam(dLamIn) {}
    ~InfoSwapMove() {}
    int i1, i2, col1, acol1, iCol1, iAcol1, col2, acol2, iCol2, iAcol2;
    double lamNow, dLam;
  };

  // Vector of (1) gluon-pair swap or (2) gluon move properties.
  vector<InfoSwapMove> infoSM;

  //----------------------------------------------------------------------

  // Initial setup: relevant gluons and coloured partons.

  inline bool setupConfig(Event& event) {

    // Reset arrays to prepare for the new event analysis.
    iGlu.clear();
    iToAllCol.resize( event.size() );
    iAllCol.clear();
    colMap.clear();
    acolMap.clear();
    infoSM.clear();
    nRec = 0;
    dLamTot = 0.;

    // Loop over all final particles. Store (fraction of) gluons.
    for (int i = 3; i < event.size(); ++i) if (event[i].isFinal()) {
      if (event[i].id() == 21 && rndmPtr->flat() < fracGluon)
        iGlu.push_back(i);

      // Store location of all colour and anticolour particles and indices.
      if (event[i].col() > 0 || event[i].acol() > 0) {
        iToAllCol[i] = iAllCol.size();
        iAllCol.push_back(i);
        if (event[i].col() > 0) colMap[event[i].col()] = i;
        if (event[i].acol() > 0) acolMap[event[i].acol()] = i;
      }
    }

    // Erase (anti)colours for (anti)junctions and skip adjacent gluons.
    for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {
      if (event.kindJunction(iJun) == 1) {
        for (int j = 0; j < 3; ++j) {
          int jCol = event.colJunction( iJun, j);
          for (map<int, int>::iterator colM = colMap.begin();
            colM != colMap.end(); ++colM)
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
          for (map<int, int>::iterator acolM = acolMap.begin();
            acolM != acolMap.end(); ++acolM)
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
    nGlu = iGlu.size();
    nCol = colMap.size();
    if (int(acolMap.size()) != nCol) {
      infoPtr->errorMsg("Error in MBReconUserHooks: map sizes do not match");
      return false;
    }
    map<int, int>::iterator colM = colMap.begin();
    map<int, int>::iterator acolM = acolMap.begin();
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
    nAllCol = iAllCol.size();
    lambdaij.resize( pow2(nAllCol) );
    int i, j;
    for (int iAC = 0; iAC < nAllCol - 1; ++iAC) {
      i = iAllCol[iAC];
      for (int jAC = iAC + 1; jAC < nAllCol; ++jAC) {
        j = iAllCol[jAC];
        lambdaij[nAllCol * iAC + jAC]
          = log(1. + m2( event[i], event[j]) / m2Ref);
      }
    }

    // Done.
    return true;

  }

  //----------------------------------------------------------------------

  // Swap gluons by lambda measure.

  inline bool doReconnectSwap(Event& event) {

    // Set up initial possible gluon swap pairs with lambda gains/losses.
    for (int iG1 = 0; iG1 < nGlu - 1; ++iG1) {
      int i1 = iGlu[iG1];
      for (int iG2 = iG1 + 1; iG2 < nGlu; ++iG2) {
        int i2 = iGlu[iG2];
        InfoSwapMove tmpSM( i1, i2);
        calcLamSwap( tmpSM, event);
        infoSM.push_back( tmpSM);
      }
    }
    int nPair = infoSM.size();

    // Keep on looping over swaps until no further negative dLambda.
    for ( int iSwap = 0; iSwap < nGlu; ++iSwap) {
      int    iPairMin = -1;
      double dLamMin  = 1e4;

      // Find lowest dLambda.
      for (int iPair = 0; iPair < nPair; ++iPair)
      if (infoSM[iPair].dLam < dLamMin) {
        iPairMin = iPair;
        dLamMin  = infoSM[iPair].dLam;
      }

      // Break if no shift below upper limit found.
      if (dLamMin > -dLamCut) break;
      ++nRec;
      dLamTot += dLamMin;
      int i1min = infoSM[iPairMin].i1;
      int i2min = infoSM[iPairMin].i2;

      // Swap the colours in the event record.
      int col1  = event[i1min].col();
      int acol1 = event[i1min].acol();
      int col2  = event[i2min].col();
      int acol2 = event[i2min].acol();
      event[i1min].cols( col2, acol2);
      event[i2min].cols( col1, acol1);

      // Swap the indices in the colour maps.
      colMap[col1]   = i2min;
      acolMap[acol1] = i2min;
      colMap[col2]   = i1min;
      acolMap[acol2] = i1min;

      // Remove already swapped pair from further consideration.
      infoSM[iPairMin] = infoSM.back();
      infoSM.pop_back();
      --nPair;

      // Update all pairs that have been affected.
      for (int iPair = 0; iPair < nPair; ++iPair) {
        InfoSwapMove& tmpSM = infoSM[iPair];
        if ( tmpSM.i1     == i1min || tmpSM.i1     == i2min
          || tmpSM.i2     == i1min || tmpSM.i2     == i2min
          || tmpSM.iCol1  == i1min || tmpSM.iCol1  == i2min
          || tmpSM.iAcol1 == i1min || tmpSM.iAcol1 == i2min
          || tmpSM.iCol2  == i1min || tmpSM.iCol2  == i2min
          || tmpSM.iAcol2 == i1min || tmpSM.iAcol2 == i2min)
          calcLamSwap( tmpSM, event);
      }
    }

    // Done.
    return true;

  }

  //----------------------------------------------------------------------

  // Calculate pair swap properties.

  inline void calcLamSwap( InfoSwapMove& tmpSM, Event& event) {

    // Colour line tracing to neighbours.
    tmpSM.col1   = event[tmpSM.i1].col();
    tmpSM.acol1  = event[tmpSM.i1].acol();
    tmpSM.iCol1  = acolMap[tmpSM.col1];
    tmpSM.iAcol1 = colMap[tmpSM.acol1];
    tmpSM.col2   = event[tmpSM.i2].col();
    tmpSM.acol2  = event[tmpSM.i2].acol();
    tmpSM.iCol2  = acolMap[tmpSM.col2];
    tmpSM.iAcol2 = colMap[tmpSM.acol2];

    // Lambda swap properties.
    double lam1c = lambda12( tmpSM.i1, tmpSM.iCol1);
    double lam1a = lambda12( tmpSM.i1, tmpSM.iAcol1);
    double lam2c = lambda12( tmpSM.i2, tmpSM.iCol2);
    double lam2a = lambda12( tmpSM.i2, tmpSM.iAcol2);
    double lam3c = lambda12( tmpSM.i1, tmpSM.iCol2);
    double lam3a = lambda12( tmpSM.i1, tmpSM.iAcol2);
    double lam4c = lambda12( tmpSM.i2, tmpSM.iCol1);
    double lam4a = lambda12( tmpSM.i2, tmpSM.iAcol1);
    if (tmpSM.col1 == tmpSM.acol2 && tmpSM.acol1 == tmpSM.col2)
       tmpSM.dLam = 2e4;
    else if (tmpSM.col1 == tmpSM.acol2)
       tmpSM.dLam = (lam3c + lam4a) - (lam1a + lam2c);
    else if (tmpSM.acol1 == tmpSM.col2)
       tmpSM.dLam = (lam3a + lam4c) - (lam1c + lam2a);
    else tmpSM.dLam = (lam3c + lam3a + lam4c + lam4a)
                   - (lam1c + lam1a + lam2c + lam2a);

  // Done.
  }

  //----------------------------------------------------------------------

  // Move gluons by lambda measure.

  inline bool doReconnectMove(Event& event) {

    // Temporary variables.
    int    iNow, colNow, acolNow, iColNow, iAcolNow, col2Now;
    double lamNow;

    // Set up initial possible gluon moves with lambda gains/losses.
    for (int iG = 0; iG < nGlu; ++iG) {

      // Gluon and its neighbours.
      iNow     = iGlu[iG];
      colNow   = event[iNow].col();
      acolNow  = event[iNow].acol();
      iColNow  = acolMap[colNow];
      iAcolNow = colMap[acolNow];

      // Addition to Lambda of gluon in current position.
      lamNow   = lambda123( iNow, iColNow, iAcolNow);

      // Loop over all colour lines where gluon could be inserted.
      for (map<int, int>::iterator colM = colMap.begin();
      colM != colMap.end(); ++colM) {
        col2Now = colM->first;

        // New container for gluon and colour line information.
        InfoSwapMove tmpSM( iNow);
        tmpSM.col1   = colNow;
        tmpSM.acol1  = acolNow;
        tmpSM.iCol1  = iColNow;
        tmpSM.iAcol1 = iAcolNow;
        tmpSM.lamNow = lamNow;
        tmpSM.col2   = col2Now;
        tmpSM.iCol2  = colMap[col2Now];
        tmpSM.iAcol2 = acolMap[col2Now];

        // Addition to Lambda if gluon inserted on line.
        tmpSM.dLam = (tmpSM.iCol2 == tmpSM.i1 || tmpSM.iAcol2 == tmpSM.i1
          || tmpSM.iCol1 == tmpSM.iAcol1) ? 2e4
          : lambda123( tmpSM.i1, tmpSM.iCol2, tmpSM.iAcol2) - tmpSM.lamNow;
        infoSM.push_back( tmpSM);
      }
    }
    int nPair = infoSM.size();

    // Keep on looping over moves until no further negative dLambda.
    for ( int iMove = 0; iMove < nGlu; ++iMove) {
      int    iPairMin = -1;
      double dLamMin  = 1e4;

      // Find lowest dLambda.
      for (int iPair = 0; iPair < nPair; ++iPair)
      if (infoSM[iPair].dLam < dLamMin) {
        iPairMin = iPair;
        dLamMin  = infoSM[iPair].dLam;
      }

      // Break if no shift below upper limit found.
      if (dLamMin > -dLamCut) break;
      ++nRec;
      dLamTot += dLamMin;

      // Partons and colours involved in move.
      InfoSwapMove& selSM = infoSM[iPairMin];
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
        InfoSwapMove& tmpSM = infoSM[iPair];
        if (tmpSM.i1 != i1Now) {
          i1Now = tmpSM.i1;
          doUpdate = false;
          if (i1Now == i1Sel || i1Now == iCol1Sel || i1Now == iAcol1Sel
            || i1Now == iCol2Sel || i1Now == iAcol2Sel) {
            colNow   = event[i1Now].col();
            acolNow  = event[i1Now].acol();
            iColNow  = acolMap[colNow];
            iAcolNow = colMap[acolNow];
            lamNow   = lambda123( i1Now, iColNow, iAcolNow);
            doUpdate = true;
          }
        }
        if (doUpdate) {
          tmpSM.col1   = colNow;
          tmpSM.acol1  = acolNow;
          tmpSM.iCol1  = iColNow;
          tmpSM.iAcol1 = iAcolNow;
          tmpSM.lamNow = lamNow;
        }
      }

      // Update info on dLambda for affected particles and colour lines.
      for (int iPair = 0; iPair < nPair; ++iPair) {
        InfoSwapMove& tmpSM = infoSM[iPair];
        int iMod = -1;
        for (int i = 0; i < 3; ++i) if (tmpSM.col2 == col2Mod[i]) iMod = i;
        if (iMod > -1) tmpSM.iCol2 = iCol2Mod[iMod];
        if (tmpSM.i1 == i1Sel || tmpSM.i1 == iCol1Sel || tmpSM.i1 == iAcol1Sel
          || tmpSM.i1 == iCol2Sel || tmpSM.i1 == iAcol2Sel || iMod > -1)
          tmpSM.dLam = (tmpSM.iCol2 == tmpSM.i1 || tmpSM.iAcol2 == tmpSM.i1
            || tmpSM.iCol1 == tmpSM.iAcol1) ? 2e4
            : lambda123( tmpSM.i1, tmpSM.iCol2, tmpSM.iAcol2) - tmpSM.lamNow;
      }

    // End of loop over gluon shifting.
    }

    // Done.
    return true;

  }

  //----------------------------------------------------------------------

  // Flip colour chains by lambda measure.

  inline bool doReconnectFlip(Event& event) {

    // Array with colour lines, and where each line begins and ends.
    vector<int> iTmp, iVec, iBeg, iEnd;

    // Grab all colour ends.
    for (int i = 3; i < event.size(); ++i)
    if (event[i].isFinal() && event[i].col() > 0 && event[i].acol() == 0) {
      iTmp.clear();
      iTmp.push_back( i);

      // Step through colour neighbours to catch system.
      int iNow = i;
      map<int, int>::iterator acolM = acolMap.find( event[iNow].col() );
      bool foundEnd = false;
      while (acolM != acolMap.end()) {
        iNow = acolM->second;
        iTmp.push_back( iNow);
        if (event[iNow].col() == 0) {
          foundEnd = true;
          break;
        }
        acolM = acolMap.find( event[iNow].col() );
      }

      // Store acceptable system, optionally including junction legs.
      if (foundEnd || flip == 2) {
        iBeg.push_back( iVec.size());
        for (int j = 0; j < int(iTmp.size()); ++j) iVec.push_back( iTmp[j]);
        iEnd.push_back( iVec.size());
      }
    }

    // Optionally search for antijunction legs: grab all anticolour ends.
    if (flip == 2) for (int i = 3; i < event.size(); ++i)
    if (event[i].isFinal() && event[i].acol() > 0 && event[i].col() == 0) {
      iTmp.clear();
      iTmp.push_back( i);

      // Step through anticolour neighbours to catch system.
      int iNow = i;
      map<int, int>::iterator colM = colMap.find( event[iNow].acol() );
      bool foundEnd = false;
      while (colM != colMap.end()) {
        iNow = colM->second;
        iTmp.push_back( iNow);
        if (event[iNow].acol() == 0) {
          foundEnd = true;
          break;
        }
        colM = colMap.find( event[iNow].acol() );
      }

      // Store acceptable system, but do not doublecount q - (n g) - qbar.
      if (!foundEnd) {
        iBeg.push_back( iVec.size());
        for (int j = 0; j < int(iTmp.size()); ++j) iVec.push_back( iTmp[j]);
        iEnd.push_back( iVec.size());
      }
    }


    // Variables for minimum search.
    int nSys = iBeg.size();
    int i1c, i1a, i2c, i2a, i1cMin, i1aMin, i2cMin, i2aMin, iSMin;
    double dLam, dLamMin;
    vector<InfoSwapMove> flipMin;

    // Loop through all system pairs.
    for (int iSys1 = 0; iSys1 < nSys - 1; ++iSys1) if (iBeg[iSys1] >= 0)
    for (int iSys2 = iSys1 + 1; iSys2 < nSys; ++iSys2) if (iBeg[iSys2] >= 0) {
      i1cMin = 0;
      i1aMin = 0;
      i2cMin = 0;
      i2aMin = 0;
      dLamMin = 1e4;

      // Loop through all possible flip locations for a pair.
      for (int j1 = iBeg[iSys1]; j1 < iEnd[iSys1] - 1; ++j1)
      for (int j2 = iBeg[iSys2]; j2 < iEnd[iSys2] - 1; ++j2) {
        i1c = iVec[j1];
        i1a = iVec[j1 + 1];
        i2c = iVec[j2];
        i2a = iVec[j2 + 1];
        dLam = lambda12( i1c, i2a) + lambda12( i2c, i1a)
             - lambda12( i1c, i1a) - lambda12( i2c, i2a);
        if (dLam < dLamMin) {
          i1cMin = i1c;
          i1aMin = i1a;
          i2cMin = i2c;
          i2aMin = i2a;
          dLamMin = dLam;
        }
      }

      // Store possible flips if low enough dLamMin.
      if (dLamMin < -dLamCut) flipMin.push_back( InfoSwapMove(
        iSys1, iSys2, i1cMin, i1aMin, i2cMin, i2aMin, dLamMin) );
    }
    int flipSize = flipMin.size();

    // Search for lowest possible flip among unused systems.
    for (int iFlip = 0; iFlip < min( nSys / 2, flipSize); ++iFlip) {
      iSMin = -1;
      dLamMin  = 1e4;
      for (int iSys12 = 0; iSys12 < flipSize; ++iSys12)
      if (flipMin[iSys12].i1 >= 0 && flipMin[iSys12].dLam < dLamMin) {
        iSMin   = iSys12;
        dLamMin = flipMin[iSys12].dLam;
      }

      // Do flip. Mark flipped systems.
      if (iSMin >= 0) {
        InfoSwapMove& flipNow = flipMin[iSMin];
        int iS1 = flipNow.i1;
        int iS2 = flipNow.i2;
        event[ flipNow.iAcol1 ].acol( event[flipNow.iCol2].col() );
        event[ flipNow.iAcol2 ].acol( event[flipNow.iCol1].col() );
        ++nRec;
        dLamTot += flipNow.dLam;
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

};

//==========================================================================


// Class for colour reconnection models specifically aimed at top decays.

class TopReconUserHooks : public UserHooks {

public:

  // Constructor and destructor.
  // mode = 0: no reconnection of tops (dummy option, does nothing);
  //      = 1: reconnect with random background gluon;
  //      = 2: reconnect with nearest (smallest-mass) background gluon;
  //      = 3: reconnect with furthest (largest-mass) background gluon;
  //      = 4: reconnect with smallest (with sign) lambda measure shift;
  //      = 5: reconnect only if reduced lamda, and then to most reduction.
  // strength: fraction of top gluons that is to be colour reconnected.
  // nList: list first nList parton classifications.
  // pTolerance: acceptable total momentum error (in reconstruction check).
  // m2Ref: squared reference mass scale for lambda measure calculation.
  // Possible variants for the future: swap with nearest in angle, not mass,
  // and/or only allow a background gluon to swap colours once.

  TopReconUserHooks(int modeIn = 0, double strengthIn = 1.) : mode(modeIn),
    strength(strengthIn) { iList = 0; nList = 0; pTolerance = 0.01;
    m2Ref = 1.;}
  ~TopReconUserHooks() {}

  // Allow colour reconnection after resonance decays (early or late)...
  virtual bool canReconnectResonanceSystems() {return true;}

  // ...which gives access to the event, for modification.
  virtual bool doReconnectResonanceSystems( int, Event& event) {

    // Return without action for relevant mode numbers.
    if (mode <= 0 || mode > 5) return true;

    // Classify coloured final partons.
    classifyFinalPartons(event);

    // Check that classification worked as expected.
    if (!checkClassification(event)) return false;

    // List first few classifications, along with the event.
    if (iList++ < nList) {
      listClassification();
      event.list();
    }

    // Perform reconnection for t and tbar in random order.
    bool tqrkFirst = (rndmPtr->flat() < 0.5);
    doReconnect( tqrkFirst, event);
    doReconnect(!tqrkFirst, event);

    // Done.
    return true;
  }

  // Return number of reconnections for current event.
  //int numberReconnections() {return nRec;}

private:

  // Mode. Counters for how many events to list. Allowed momentum error.
  int    mode, iList, nList, nRec;
  double strength, pTolerance, m2Ref;

  // Arrays of (indices of) final partons from different sources.
  // So far geared towards t -> b W decays only.
  vector<int> iBqrk, iWpos, iTqrk, iBbar, iWneg, iTbar, iRest;

  // Maps telling where all colours and anticolours are stored.
  map<int, int> colMap, acolMap;

  //----------------------------------------------------------------------

  // Classify all coloured partons at the end of showers by origin.
  // Note: for now only t -> b W is fully classified.

  inline bool classifyFinalPartons(Event& event) {

    // Reset arrays to prepare for the new event analysis.
    iBqrk.clear();
    iWpos.clear();
    iTqrk.clear();
    iBbar.clear();
    iWneg.clear();
    iTbar.clear();
    iRest.clear();
    colMap.clear();
    acolMap.clear();
    nRec = 0;

    // Loop over all final particles. Tag coloured ones.
    for (int i = 3; i < event.size(); ++i) if (event[i].isFinal()) {
      bool hasCol = (event[i].colType() != 0);

      // Set up to find where each parton comes from.
      bool fsrFromT = false;
      bool fromTqrk = false;
      bool fromBqrk = false;
      bool fromWpos = false;
      bool fromTbar = false;
      bool fromBbar = false;
      bool fromWneg = false;

      // Identify current particle.
      int iNow  = i;
      int idOld = 0;
      do {
        int idNow = event[iNow].id();

        // Exclude FSR of gluons/photons from the t quark proper.
        if (abs(idNow) == 6 && (idOld == 21 || idOld == 22)) fsrFromT = true;

        // Check if current particle matches any of the categories.
        else if (idNow ==   6) fromTqrk = true;
        else if (idNow ==   5) fromBqrk = true;
        else if (idNow ==  24) fromWpos = true;
        else if (idNow ==  -6) fromTbar = true;
        else if (idNow ==  -5) fromBbar = true;
        else if (idNow == -24) fromWneg = true;

        // Step up through the history to the very top.
        iNow  = event[iNow].mother1();
        idOld = idNow;
      } while (iNow > 2 && !fsrFromT);

      // Bookkeep where the parton comes from. Note that b quarks also
      // can appear in W decays, so order of checks is relevant.
      if      (fromTqrk && fromWpos && hasCol) iWpos.push_back(i);
      else if (fromTqrk && fromBqrk && hasCol) iBqrk.push_back(i);
      else if (fromTqrk)                       iTqrk.push_back(i);
      else if (fromTbar && fromWneg && hasCol) iWneg.push_back(i);
      else if (fromTbar && fromBbar && hasCol) iBbar.push_back(i);
      else if (fromTbar)                       iTbar.push_back(i);
      else if (hasCol)                         iRest.push_back(i);

      // Store location of all colour and anticolour indices.
      if (hasCol && (mode == 4 || mode == 5)) {
        if (event[i].col() > 0) colMap[event[i].col()] = i;
        if (event[i].acol() > 0) acolMap[event[i].acol()] = i;
      }
    }

    // So far method always returns true.
    return true;

  }

  //----------------------------------------------------------------------

  // Check that classification worked by summing up partons to t/tbar.

  inline bool checkClassification(Event& event) {

    // Find final copy of t and tbar quarks.
    int iTqrkLoc = 0;
    int iTbarLoc = 0;
    for (int i = 3; i < event.size(); ++i) {
      if(event[i].id() ==  6) iTqrkLoc = i;
      if(event[i].id() == -6) iTbarLoc = i;
    }

    // Four-momentum of t minus all its decay products.
    Vec4 tqrkDiff = event[iTqrkLoc].p();
    for (int i = 0; i < int(iBqrk.size()); ++i)
      tqrkDiff -= event[iBqrk[i]].p();
    for (int i = 0; i < int(iWpos.size()); ++i)
      tqrkDiff -= event[iWpos[i]].p();
    for (int i = 0; i < int(iTqrk.size()); ++i)
      tqrkDiff -= event[iTqrk[i]].p();

    // Four-momentum of tbar minus all its decay products.
    Vec4 tbarDiff = event[iTbarLoc].p();
    for (int i = 0; i < int(iBbar.size()); ++i)
      tbarDiff -= event[iBbar[i]].p();
    for (int i = 0; i < int(iWneg.size()); ++i)
      tbarDiff -= event[iWneg[i]].p();
    for (int i = 0; i < int(iTbar.size()); ++i)
      tbarDiff -= event[iTbar[i]].p();

    // Print difference vectors and event if sum deviation is too big.
    double totErr = abs(tqrkDiff.px()) + abs(tqrkDiff.py())
      + abs(tqrkDiff.pz()) + abs(tqrkDiff.e()) + abs(tbarDiff.px())
      + abs(tbarDiff.py()) + abs(tbarDiff.pz()) + abs(tqrkDiff.e());
    if (totErr > pTolerance) {
      infoPtr->errorMsg("Error in TopReconUserHooks::checkClassification");
      cout << "\n Error in t/tbar daughter search: \n t    difference "
           << tqrkDiff << " tbar difference "<< tbarDiff;
      listClassification();
      event.list();
    }

    // Done.
    return (totErr < pTolerance);
  }

  //----------------------------------------------------------------------

  // Print how final-state (mainly coloured) particles were classified.

  inline void listClassification() {

    cout << "\n Final-state coloured partons classified by source: ";
    cout << "\n From Bqrk:";
    for (int i = 0; i < int(iBqrk.size()); ++i) cout << "  " << iBqrk[i];
    cout << "\n From Wpos:";
    for (int i = 0; i < int(iWpos.size()); ++i) cout << "  " << iWpos[i];
    cout << "\n From Tqrk:";
    for (int i = 0; i < int(iTqrk.size()); ++i) cout << "  " << iTqrk[i];
    cout << "\n From Bbar:";
    for (int i = 0; i < int(iBbar.size()); ++i) cout << "  " << iBbar[i];
    cout << "\n From Wneg:";
    for (int i = 0; i < int(iWneg.size()); ++i) cout << "  " << iWneg[i];
    cout << "\n From Tbar:";
    for (int i = 0; i < int(iTbar.size()); ++i) cout << "  " << iTbar[i];
    cout << "\n From Rest:";
    for (int i = 0; i < int(iRest.size()); ++i) {
      cout << "  " << iRest[i];
      if (i%20 == 19 && i + 1 != int(iRest.size())) cout << "\n           ";
    }
    cout << endl;
  }

  //----------------------------------------------------------------------

  // Reconnect gluons either from t or from tbar quark.

  inline bool doReconnect(bool doTqrk, Event& event) {

    // Gather coloured decay products either of t or of tbar.
    vector<int> iTdec;
    if (doTqrk) {
      for (int i = 0; i < int(iBqrk.size()); ++i) iTdec.push_back(iBqrk[i]);
      for (int i = 0; i < int(iWpos.size()); ++i) iTdec.push_back(iWpos[i]);
    } else {
      for (int i = 0; i < int(iBbar.size()); ++i) iTdec.push_back(iBbar[i]);
      for (int i = 0; i < int(iWneg.size()); ++i) iTdec.push_back(iWneg[i]);
    }

    // Extract the gluons from the t quark decay.
    vector<int> iGT;
    for (int i = 0; i < int(iTdec.size()); ++i) {
      int colNow  = event[iTdec[i]].col();
      int acolNow = event[iTdec[i]].acol();
      if (colNow > 0 && acolNow > 0) iGT.push_back(iTdec[i]);
    }
    int nGT    = iGT.size();

    // Randomize their stored order.
    if (nGT > 1) for (int i = 0; i < nGT; ++i) {
      int j = min( int(nGT * rndmPtr->flat()), nGT - 1 );
      swap( iGT[i], iGT[j]);
    }

    // Also extract the rest of the gluons in the event.
    vector<int> iGR;
    for (int i = 0; i < int(iRest.size()); ++i) {
      int colNow = event[iRest[i]].col();
      int acolNow = event[iRest[i]].acol();
      if (colNow > 0 && acolNow > 0) iGR.push_back(iRest[i]);
    }
    int nGR    = iGR.size();
    int iR, colT, acolT, iColT, iAcolT, colR, acolR, iColR, iAcolR;
    double mTR2, mTR2now, dLam, lamT, lamNow, lamRec;

    // Loop through all top gluons; study fraction given by strength.
    if (nGT > 0 && nGR > 0)
    for (int iT = 0; iT < nGT; ++iT) {
      if (strength < rndmPtr->flat()) continue;

      // Pick random gluon from rest of event.
      if (mode == 1) iR = min( int(nGR * rndmPtr->flat()), nGR - 1 );

      // Find gluon from rest with lowest or highest invariant mass.
      else if (mode < 4) {
        iR   = 0;
        mTR2 = m2( event[iGT[iT]], event[iGR[iR]]);
        for (int ii = 1; ii < nGR; ++ii) {
          mTR2now = m2( event[iGT[iT]], event[iGR[ii]]);
          if (mode == 2 && mTR2now < mTR2) {iR = ii; mTR2 = mTR2now;}
          if (mode == 3 && mTR2now > mTR2) {iR = ii; mTR2 = mTR2now;}
        }

      // Find gluon from rest with smallest lambda value shift.
      } else {
        iR     = -1;
        dLam   = 1e10;
        colT   = event[iGT[iT]].col();
        acolT  = event[iGT[iT]].acol();
        iColT  = acolMap[colT];
        iAcolT = colMap[acolT];
        lamT   = log(1. + m2( event[iGT[iT]], event[iColT]) / m2Ref)
               + log(1. + m2( event[iGT[iT]], event[iAcolT]) / m2Ref);
        for (int ii = 0; ii < nGR; ++ii) {
          colR   = event[iGR[ii]].col();
          acolR  = event[iGR[ii]].acol();
          iColR  = acolMap[colR];
          iAcolR = colMap[acolR];
          lamNow = lamT
                 + log(1. + m2( event[iGR[ii]], event[iColR]) / m2Ref)
                 + log(1. + m2( event[iGR[ii]], event[iAcolR]) / m2Ref);
          lamRec = log(1. + m2( event[iGT[iT]], event[iColR]) / m2Ref)
                 + log(1. + m2( event[iGT[iT]], event[iAcolR]) / m2Ref)
                 + log(1. + m2( event[iGR[ii]], event[iColT]) / m2Ref)
                 + log(1. + m2( event[iGR[ii]], event[iAcolT]) / m2Ref);
          if (lamRec - lamNow < dLam) {iR = ii; dLam = lamRec - lamNow;}
        }
      }
      if (mode == 5 && dLam > 0.) continue;

      // Swap top and rest gluon colour and anticolour.
      ++nRec;
      swapCols( iGT[iT], iGR[iR], event);
    }

    // Done.
    return true;

  }

  //----------------------------------------------------------------------

  // Swap colours and/or anticolours in the event listing.

  inline void swapCols( int i, int j, Event& event) {

    // Swap the colours in the event record.
    int coli  = event[i].col();
    int acoli = event[i].acol();
    int colj  = event[j].col();
    int acolj = event[j].acol();
    event[i].cols( colj, acolj);
    event[j].cols( coli, acoli);

    // Swap the indices in the colour maps.
    colMap[coli]   = j;
    acolMap[acoli] = j;
    colMap[colj]   = i;
    acolMap[acolj] = i;

  }

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_ColourReconnectionHooks_H
