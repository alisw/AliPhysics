// ColourTracing.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// ColourReconnection class.

// Setup the list of colours, this is needed later for finding colour chains.

#include "Pythia8/ColourTracing.h"

namespace Pythia8 {

//==========================================================================

// The ColourTracing class.

//--------------------------------------------------------------------------

// Find all final coloured and anticoloured partons.

bool ColourTracing::setupColList(Event& event) {

  iColEnd.resize(0);
  iAcolEnd.resize(0);
  iColAndAcol.resize(0);
  for (int i = 0; i < event.size(); ++i)
    if (event[i].isFinal()) {
      if (event[i].col() > 0 && event[i].acol() > 0) iColAndAcol.push_back(i);
      else if (event[i].col() > 0) iColEnd.push_back(i);
      else if (event[i].acol() > 0) iAcolEnd.push_back(i);
      // Colour sextets have additional tags (store with negative numbers).
      if (event[i].col() < 0) iAcolEnd.push_back(-i);
      else if (event[i].acol() < 0) iColEnd.push_back(-i);
    }

  // Return true if zero particles were found.
  if (int(iColEnd.size()) == 0 && int(iAcolEnd.size()) == 0 &&
      int(iColAndAcol.size()) == 0) return true;
  else return false;

}

//--------------------------------------------------------------------------

// Trace a colour line, from an anticolour to a colour.

bool ColourTracing::traceFromAcol(int indxCol, Event& event, int iJun,
  int iCol, vector<int>& iParton) {

  // Junction kind, if any.
  int kindJun = (iJun >= 0) ? event.kindJunction(iJun) : 0;

  // Begin to look for a matching colour.
  int loop = 0;
  int loopMax = iColAndAcol.size() + 2;
  bool hasFound = false;
  do {
    ++loop;
    hasFound= false;

    // First check list of matching colour ends.
    // Also check for sextets (negative anticolour tag = extra colour tag).
    for (int i = 0; i < int(iColEnd.size()); ++i) {
      if (event[ abs(iColEnd[i]) ].col() == indxCol
          || event[ abs(iColEnd[i]) ].acol() == -indxCol) {
        iParton.push_back( abs(iColEnd[i]) );
        indxCol = 0;
        iColEnd[i] = iColEnd.back();
        iColEnd.pop_back();
        hasFound = true;
        break;
      }
    }

    // Then check list of intermediate gluons.
    if (!hasFound)
    for (int i = 0; i < int(iColAndAcol.size()); ++i)
    if (event[ iColAndAcol[i] ].col() == indxCol) {
      iParton.push_back( iColAndAcol[i] );
      // Update to new colour. Remove gluon.
      indxCol = event[ iColAndAcol[i] ].acol();
      if (kindJun > 0) event.endColJunction(iJun, iCol, indxCol);
      iColAndAcol[i] = iColAndAcol.back();
      iColAndAcol.pop_back();
      hasFound = true;
      break;
    }

    // Check opposite-sign junction colours.
    if (!hasFound)
    for (int iAntiJun = 0; iAntiJun < event.sizeJunction(); ++iAntiJun)
      if (iAntiJun != iJun && event.kindJunction(iAntiJun) % 2 == 0)
        for (int iColAnti = 0; iColAnti < 3; ++iColAnti)
          if (event.colJunction(iAntiJun, iColAnti) == indxCol) {
            iParton.push_back( -(10 + 10 * iAntiJun + iColAnti) );
            indxCol = 0;
            hasFound = true;
            break;
          }

    // In a pinch, check list of opposite-sign junction end colours.
    // Store in iParton list as -(10 + 10 * iAntiJun + iLeg).
    // This is for J-g-...-g-J connections; where instead of running both ways,
    // the second time we just store the two junctions.
    if (!hasFound && kindJun % 2 == 1 && event.sizeJunction() > 1)
    for (int iAntiJun = 0; iAntiJun < event.sizeJunction(); ++iAntiJun)
      if (iAntiJun != iJun && event.kindJunction(iAntiJun) % 2 == 0)
        for (int iColAnti = 0; iColAnti < 3; ++iColAnti)
          if (event.endColJunction(iAntiJun, iColAnti) == indxCol) {
            iParton.push_back( -(10 + 10 * iAntiJun + iColAnti) );
            indxCol = 0;
            hasFound = true;
            break;
          }

    // Keep on tracing via gluons until reached end of leg.
  } while (hasFound && indxCol > 0 && loop < loopMax);

  // Something went wrong in colour tracing.
  if (!hasFound || loop == loopMax) {
    infoPtr->errorMsg("Error in ColourTracing::traceFromAcol: "
      "colour tracing failed");
    return false;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Trace a colour line, from a colour to an anticolour.

bool ColourTracing::traceFromCol(int indxCol, Event& event, int iJun,
  int iCol, vector<int>& iParton) {

  // If none specified, select next colour tag from back of list.
  if (iJun < 0  && iCol < 0) {
    int iColEndBack = iColEnd.back();
    if (iColEndBack > 0) indxCol = event[iColEnd.back()].col();
    // Negative index implies extra (sextet) colour tag in anticolour slot.
    else                 indxCol = -event[-iColEnd.back()].acol();
    iParton.push_back(iColEnd.back());
    iColEnd.pop_back();
  }

  // Junction kind, if any.
  int kindJun = (iJun >= 0) ? event.kindJunction(iJun) : 0;

  // Begin to look for a matching anticolour.
  int loop = 0;
  int loopMax = iColAndAcol.size() + 2;
  bool hasFound = false;
  do {
    ++loop;
    hasFound= false;

    // First check list of matching anticolour ends.
    // Also check for sextets (negative colour tag = extra anticolour tag).
    for (int i = 0; i < int(iAcolEnd.size()); ++i) {
      if (event[ abs(iAcolEnd[i]) ].acol() == indxCol
          || event[ abs(iAcolEnd[i]) ].col() == -indxCol) {
        iParton.push_back( abs(iAcolEnd[i]) );
        indxCol = 0;
        iAcolEnd[i] = iAcolEnd.back();
        iAcolEnd.pop_back();
        hasFound = true;
        break;
      }
    }

    // Then check list of intermediate gluons.
    if (!hasFound)
    for (int i = 0; i < int(iColAndAcol.size()); ++i)
    if (event[ iColAndAcol[i] ].acol() == indxCol) {
      iParton.push_back( iColAndAcol[i] );

      // Update to new colour. Remove gluon.
      indxCol = event[ iColAndAcol[i] ].col();
      if (kindJun > 0) event.endColJunction(iJun, iCol, indxCol);
      iColAndAcol[i] = iColAndAcol.back();
      iColAndAcol.pop_back();
      hasFound = true;
      break;
    }

    // Check opposite-sign junction colours.
    if (!hasFound)
      for (int iAntiJun = 0; iAntiJun < event.sizeJunction(); ++iAntiJun)
        if (iAntiJun != iJun && event.kindJunction(iAntiJun) %2 == 1)
          for (int iColAnti = 0; iColAnti < 3; ++iColAnti)
            if (event.colJunction(iAntiJun, iColAnti) == indxCol) {
              iParton.push_back( -(10 + 10 * iAntiJun + iColAnti) );
              indxCol = 0;
              hasFound = true;
              break;
            }

    // In a pinch, check list of opposite-sign junction end colours.
    // Store in iParton list as -(10 + 10 * iAntiJun + iAntiLeg).
    // This is for J-g-...-g-J connections; where instead of running both ways,
    // the second time we just store the two junctions.
    if (!hasFound && kindJun % 2 == 0 && event.sizeJunction() > 1)
      for (int iAntiJun = 0; iAntiJun < event.sizeJunction(); ++iAntiJun)
        if (iAntiJun != iJun && event.kindJunction(iAntiJun) %2 == 1)
          for (int iColAnti = 0; iColAnti < 3; ++iColAnti)
            if (event.endColJunction(iAntiJun, iColAnti) == indxCol) {
              iParton.push_back( -(10 + 10 * iAntiJun + iColAnti) );
              indxCol = 0;
              hasFound = true;
              break;
            }

  // Keep on tracing via gluons until reached end of leg.
  } while (hasFound && indxCol > 0 && loop < loopMax);

  // Something went wrong in colour tracing.
  if (!hasFound || loop == loopMax) {
    infoPtr->errorMsg("Error in ColourTracing::traceFromCol: "
      "colour tracing failed");
    return false;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Trace a colour loop, from a colour back to the anticolour of the same.

bool ColourTracing::traceInLoop(Event& event, vector<int>& iParton) {

  // Add starting gluon.
  iParton.push_back( iColAndAcol[0] );
  int indxCol = event[ iColAndAcol[0] ].col();
  int indxAcol = event[ iColAndAcol[0] ].acol();
  iColAndAcol[0] = iColAndAcol.back();
  iColAndAcol.pop_back();

  // Move around until back where begun.
  int loop = 0;
  int loopMax = iColAndAcol.size() + 2;
  bool hasFound = false;
  do {
    ++loop;
    hasFound= false;

    // Check list of gluons.
    for (int i = 0; i < int(iColAndAcol.size()); ++i)
      if (event[ iColAndAcol[i] ].acol() == indxCol) {
        iParton.push_back( iColAndAcol[i] );
        indxCol = event[ iColAndAcol[i] ].col();
        iColAndAcol[i] = iColAndAcol.back();
        iColAndAcol.pop_back();
        hasFound = true;
        break;
    }
  } while (hasFound && indxCol != indxAcol && loop < loopMax);

  // Something went wrong in colour tracing.
  if (!hasFound || loop == loopMax) {
    infoPtr->errorMsg("Error in ColourTracing::traceInLoop: "
      "colour tracing failed");

    return false;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Get junction chains, where the junctions are directly connected.

vector<vector<int > > ColourTracing::getJunChains(Event& event) {

  // Make list of junction chains and help array.
  vector<vector<int> > junChains;
  vector<bool> usedJuncs(event.sizeJunction(),false);

  // Loop over junctions.
  for (int i = 0; i < event.sizeJunction(); ++i) {
    if (usedJuncs[i])
      continue;
    std::list<int> curJun;
    vector<int> junList;
    usedJuncs[i] = true;
    curJun.push_back(i);
    junList.push_back(i);

    // Keep looping over connected junctions until no new junctions are found.
    while (!curJun.empty()) {
      for (int iLeg = 0;iLeg < 3; ++iLeg)
        for (int j = 0;j < event.sizeJunction(); ++j) {
          if (usedJuncs[j])
            continue;
          for (int jLeg = 0;jLeg < 3; ++jLeg) {
            if (event.colJunction(curJun.front(),iLeg) ==
                event.colJunction(j,jLeg)) {
              curJun.push_back(j);
              junList.push_back(j);
              usedJuncs[j] = true;
              break;
            }
          }
        }
      curJun.pop_front();
    }
    junChains.push_back(junList);
  }
  return junChains;

}

//==========================================================================

} // end namespace Pythia8
