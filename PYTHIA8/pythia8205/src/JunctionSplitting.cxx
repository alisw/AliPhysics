// JunctionSplitting.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// JunctionSplitting class.

// Setup the list of colours, this is needed later for finding colour chains.

#include "Pythia8/JunctionSplitting.h"

namespace Pythia8 {

//==========================================================================

// The JunctionSplitting class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// For breaking J-J string, pick a Gamma by taking a step with fictitious mass.
const double JunctionSplitting::JJSTRINGM2MAX  = 25.;
const double JunctionSplitting::JJSTRINGM2FRAC = 0.1;

// Iterate junction rest frame boost until convergence or too many tries.
const double JunctionSplitting::CONVJNREST     = 1e-5;
const int JunctionSplitting::NTRYJNREST        = 20;

// Typical average transvere primary hadron mass <mThad>.
const double JunctionSplitting::MTHAD          = 0.9;

//--------------------------------------------------------------------------

// Initialize the class and all the created classes.

void JunctionSplitting::init( Info* infoPtrIn, Settings& settings,
  Rndm* rndmPtrIn, ParticleData* particleDataPtrIn) {

  infoPtr = infoPtrIn;
  rndmPtr = rndmPtrIn;

  // Initialize
  colTrace.init(infoPtrIn);
  stringLength.init(infoPtrIn, settings);

  // Initialize auxiliary fragmentation classes.
  flavSel.init(settings, rndmPtr);
  pTSel.init(settings, *particleDataPtrIn, rndmPtr);
  zSel.init(settings, *particleDataPtrIn, rndmPtr);

  // Initialize string and ministring fragmentation.
  stringFrag.init(infoPtr, settings, particleDataPtrIn, rndmPtr,
    &flavSel, &pTSel, &zSel);

  // For junction processing.
  eNormJunction     = settings.parm("StringFragmentation:eNormJunction");
  allowDoubleJunRem = settings.flag("ColourReconnection:allowDoubleJunRem");
}

//--------------------------------------------------------------------------

// Check that all colours are connected in physical way. Also split
// junction pairs, such that the hadronization can handle the configuration.

bool JunctionSplitting::checkColours( Event& event) {

  // Not really a colour check, but a check all numbers are valid.
  for (int i = 0; i < event.size(); ++i)
    if (abs(event[i].px()) >= 0. && abs(event[i].py()) >= 0.
      && abs(event[i].pz()) >= 0.  && abs(event[i].e()) >= 0.
       && abs(event[i].m()) >= 0.);
    else {
       infoPtr->errorMsg("Error in Pythia::CheckColours: "
        "not-a-number energy/momentum/mass");
       return false;
    }

  // Check if any singlet gluons were made, and if so return false.
  for (int i  = 0; i < event.size(); ++i) {
    if (event[i].isFinal() && event[i].col() != 0 &&
        event[i].col() == event[i].acol()) {
      infoPtr->errorMsg("Warning in JunctionSplitting::CheckColours:"
      "Made a gluon colour singlet. Redoing colour configuration");
      return false;
    }
  }

  // Need to try and split junction structures.
  colTrace.setupColList(event);
  vector<int> iParton;
  vector<vector <int > > iPartonJun, iPartonAntiJun;
  getPartonLists(event, iPartonJun, iPartonAntiJun);

  // Try to split up the junction chains by splitting gluons
  if (!splitJunGluons(event, iPartonJun, iPartonAntiJun) ) {
    infoPtr->errorMsg("Warning in JunctionSplitting::CheckColours:"
      "Not possible to split junctions. Making new colour configuration");
    return false;
  }

  // Remove junctions if more than 2 are connected.
  if (!splitJunChains(event) ) {
    infoPtr->errorMsg("Warning in JunctionSplitting::CheckColours:"
      "Not possible to split junctions. Making new colour configuration");
    return false;
  }

  // Split up junction pairs.
  getPartonLists(event, iPartonJun, iPartonAntiJun);
  if (!splitJunPairs(event, iPartonJun, iPartonAntiJun) ) {
    infoPtr->errorMsg("Warning in JunctionSplitting::CheckColours:"
      "Not possible to split junctions. Making new colour configuration");
    return false;
  }

  // Done checking.
  return true;
}

//--------------------------------------------------------------------------

// Split connected junction chains into separated, mainly by splitting gluons
// into q-qbar pairs. If the junctions are directly connected
// other methods are applied.

bool JunctionSplitting::splitJunGluons(Event& event,
  vector<vector< int > >& iPartonJun, vector<vector< int > >& iPartonAntiJun) {

  // Loop over all junctions and all junction legs.
  for (int iJun = 0; iJun < int(iPartonJun.size()); ++iJun) {

    // Fill in vector of the legs content.
    vector<vector <int> > iJunLegs;
    iJunLegs.resize(3);
    int leg = -1;
    for (int i = 0; i < int(iPartonJun[iJun].size()); ++i) {
      if ( iPartonJun[iJun][i]/10 == iPartonJun[iJun][0]/10) ++leg;
      iJunLegs[leg].push_back(iPartonJun[iJun][i]);
    }

    // Loop over legs.
    for (int i = 0;i < int(iJunLegs.size()); ++i) {
      // If it is not connected to another junction, no need to do anything.
      if (iJunLegs[i].back() > 0)
        continue;
      int identJun = iJunLegs[i][0];
      // If no gluons in between two junctions, not possible to do anything.
      if (iJunLegs[i].size() == 2)
        continue;

      int identAntiJun = 0, iAntiLeg = -1;

      // Pick a new quark at random; for simplicity no diquarks.
      int colQ = 0, acolQ = 0;
      int idQ = int(rndmPtr->flat() * 3) + 1;

      // If a single gluon in between the two junctions, change it to a
      // quark-anti quark system.
      if ( iJunLegs[i].size() == 3) {

        // Store the new q qbar pair, sharing gluon colour and momentum.
        colQ = event[ iJunLegs[i][1] ].col();
        acolQ = event[ iJunLegs[i][1] ].acol();
        Vec4 pQ = 0.5 * event[ iJunLegs[i][1] ].p();
        double mQ = 0.5 * event[ iJunLegs[i][1] ].m();
        int iQ = event.append( idQ, 75, iJunLegs[i][1], 0, 0, 0, colQ, 0,
          pQ, mQ );
        int iQbar = event.append( -idQ, 75, iJunLegs[i][1], 0, 0, 0, 0, acolQ,
          pQ, mQ );
        
        // Mark split gluon.
        event[ iJunLegs[i][1] ].statusNeg();
        event[ iJunLegs[i][1] ].daughters( iQ, iQbar);

        // Update junction and anti junction list.
        identAntiJun = iJunLegs[i].back();
        int iOld = iJunLegs[i][1];
        bool erasing = false;
        for (int j = 0; j < int(iPartonJun[iJun].size()); ++j) {
          if (iPartonJun[iJun][j] == iOld)
            erasing = true;
          if (iPartonJun[iJun][j] == identAntiJun) {
            iPartonJun[iJun][j] = iQ;
            break;
          }
          if (erasing) {
            iPartonJun[iJun].erase(iPartonJun[iJun].begin() + j);
            --j;
          }
        }

        // Find the connected anti junction from the list of anti junctions.
        int iAntiJun = -1;
        for (int j = 0; j < int(iPartonAntiJun.size()); j++)
          if ( iPartonAntiJun[j][0]/10 == identAntiJun/10) {
            iAntiJun = j;
            break;
          }
        // If no anti junction found, something went wrong earlier.
        if (iAntiJun == -1) {
           infoPtr->errorMsg("Warning in JunctionSplitting::SplitJunChain:"
                             "Something went wrong in finding anti junction");
           return false;
        }

        // Update the anti junction list.
        erasing = false;
        for (int j = 0; j < int(iPartonAntiJun[iAntiJun].size()); ++j) {
          if ( iPartonAntiJun[iAntiJun][j] / 10 == identAntiJun / 10)
            iAntiLeg++;
          if ( iPartonAntiJun[iAntiJun][j] == identJun) {
            iPartonAntiJun[iAntiJun][j] = iQbar;
            break;
          }
        }
      }
      // If more than a single gluon, decide depending on mass.
      else if (iJunLegs[i].size() > 3) {
        // Evaluate mass-squared for all adjacent gluon pairs.
        vector<double > m2Pair;
        double m2Sum = 0.;

        for (int j = 1; j < int(iJunLegs[i].size()) - 2; ++j) {
          double m2Now = 0.5 * event[ iJunLegs[i][j] ].p()
            * event[ iJunLegs[i][j + 1] ].p();
          m2Pair.push_back(m2Now);
          m2Sum += m2Now;
        }

        // Pick breakup region with probability proportional to mass-squared.
        double m2Reg = m2Sum * rndmPtr->flat();
        int iReg = -1;
        do m2Reg -= m2Pair[++iReg];
        while (m2Reg > 0. && iReg < int(iJunLegs[i].size()) - 1);
        m2Reg = m2Pair[iReg];

        // increase iReg with one, since it should not point towards itself.
        iReg++;

        // Pick breaking point of string in chosen region (symmetrically).
        double m2Temp = min( JJSTRINGM2MAX, JJSTRINGM2FRAC * m2Reg);
        double xPos = 0.5;
        double xNeg = 0.5;
        do {
          double zTemp = zSel.zFrag( idQ, 0, m2Temp);
          xPos = 1. - zTemp;
          xNeg = m2Temp / (zTemp * m2Reg);
        } while (xNeg > 1.);
        if (rndmPtr->flat() > 0.5) swap(xPos, xNeg);

        // Pick up two "mother" gluons of breakup. Mark them decayed.
        Particle& gJun = event[ iJunLegs[i][iReg] ];
        Particle& gAnti = event[ iJunLegs[i][iReg + 1] ];
        gJun.statusNeg();
        gAnti.statusNeg();
        int dau1 = event.size();
        gJun.daughters(dau1, dau1 + 3);
        gAnti.daughters(dau1, dau1 + 3);
        int mother1 = min( iJunLegs[i][iReg], iJunLegs[i][iReg + 1]);
        int mother2 = max( iJunLegs[i][iReg], iJunLegs[i][iReg + 1]);

        // Need to store variables, since it is not safe to use references
        // with append.
        int gJunCol   = gJun.col();
        int gJunAcol  = gJun.acol();
        int gAntiAcol = gAnti.acol();
        Vec4 gJunP    = gJun.p();
        Vec4 gAntiP   = gAnti.p();
        double gJunM  = gJun.m();
        double gAntiM = gAnti.m();

        // Can keep one of old colours but need one new so unambiguous.
        colQ          = gJunAcol;
        acolQ         = event.nextColTag();

         // Store copied gluons with reduced momenta.
        int iGjun = event.append( 21, 75, mother1, mother2, 0, 0,
          gJunCol, gJunAcol, (1. - 0.5 * xPos) * gJunP,
          (1. - 0.5 * xPos) * gJunM);
        event.append( 21, 75, mother1, mother2, 0, 0,
          acolQ, gAntiAcol, (1. - 0.5 * xNeg) * gAntiP,
          (1. - 0.5 * xNeg) * gAntiM);

        // Store the new q qbar pair with remaining momenta.
        int iQ = event.append( idQ, 75, mother1, mother2, 0, 0,
          colQ, 0, 0.5 * xNeg * gAntiP, 0.5 * xNeg * gAntiM );
        int iQbar = event.append( -idQ, 75, mother1, mother2, 0, 0,
          0, acolQ, 0.5 * xPos * gJunP, 0.5 * xPos * gJunM );

        // Update the list of junctions to reflect the splitting.
        identAntiJun = iJunLegs[i].back();
        bool erasing = false;
        for (int j = 0; j < int(iPartonJun[iJun].size()); ++j) {
          if (iPartonJun[iJun][j] == mother1 ||
              iPartonJun[iJun][j] == mother2)
            erasing = true;
        
          if ( iPartonJun[iJun][j] == identAntiJun) {
            iPartonJun[iJun][j] = iQ;
            iPartonJun[iJun].insert(iPartonJun[iJun].begin() + j, iGjun);
            break;
          }
          if (erasing) {
            iPartonJun[iJun].erase(iPartonJun[iJun].begin() + j);
            j--;
          }
        }

        // Find the connected anti junction from the list of anti junctions.
        int iAntiJun = -1;
        for (int j = 0; j < int(iPartonAntiJun.size());j++)
          if ( iPartonAntiJun[j][0]/10 == identAntiJun/10) {
            iAntiJun = j;
            break;
          }
        // If no anti junction found, something went wrong earlier.
        if (iAntiJun == -1) {
           infoPtr->errorMsg("Error in JunctionSplitting::SplitJunChain:"
                             "Something went wrong in finding anti junction");
           return false;
        }

        // Update the anti junction list to reflect the splitting
        for (int j = 0; j < int(iPartonAntiJun[iAntiJun].size()); ++j) {
          if ( iPartonAntiJun[iAntiJun][j] / 10 == identAntiJun / 10)
            iAntiLeg++;
          if (iPartonAntiJun[iAntiJun][j] == identJun) {
            iPartonAntiJun[iAntiJun][j] = iQbar;
            break;
          }
        }       
      }

      // Update end colours for both g -> q qbar and g g -> g g q qbar.
      event.endColJunction((-identJun)/10 - 1, i, colQ);
      event.endColJunction((-identAntiJun)/10 - 1, iAntiLeg, acolQ);
    }
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Fix chains that contain more than two junctions.
// This is done by removing the minimum needed amount of junctions.
// Might need to make choice based on String length, now randomly chosen.

bool JunctionSplitting::splitJunChains(Event& event) {

  // Get junction chains.
  event.saveJunctionSize();
  vector<vector<int> > junChains = colTrace.getJunChains(event);

  // Remove junctions.
  vector<int> junRem;
  for (int i = 0; i < int(junChains.size()); ++i) {
    if (junChains[i].size() < 3)
      continue;

    vector<int> cols, acols;
    for (int j = 0; j < int(junChains[i].size()); ++j) {

      junRem.push_back(junChains[i][j]);
      if (event.kindJunction(junChains[i][j]) % 2 == 0)
        for (int jLeg = 0; jLeg < 3; ++jLeg)
          acols.push_back(event.colJunction(junChains[i][j],jLeg));
      else
        for (int jLeg = 0; jLeg < 3; ++jLeg)
          cols.push_back(event.colJunction(junChains[i][j],jLeg));
    }

    for (int j = 0; j < int(cols.size()); ++j)
      for (int k = 0; k < int(acols.size()); ++k)
        if (cols[j] == acols[k]) {
          cols.erase(cols.begin() + j);
          acols.erase(acols.begin() + k);
          j--;
          break;
        }

    // Find junctions if we have more colours than anti colours
    while (cols.size() > acols.size()) {
      int i1 = int(rndmPtr->flat() *cols.size());
      int col1 = cols[i1];
      cols.erase(cols.begin() + i1);
      int i2 = int(rndmPtr->flat() *cols.size());
      int col2 = cols[i2];
      cols.erase(cols.begin() + i2);
      int i3 = int(rndmPtr->flat() *cols.size());
      int col3 = cols[i3];
      cols.erase(cols.begin() + i3);
      event.appendJunction(1, col1, col2, col3);
    }

    // Find junctions if we have more colours than anti colours
    while (acols.size() > cols.size()) {
      int i1 = int(rndmPtr->flat() *acols.size());
      int acol1 = acols[i1];
      acols.erase(acols.begin() + i1);
      int i2 = int(rndmPtr->flat() *acols.size());
      int acol2 = acols[i2];
      acols.erase(acols.begin() + i2);
      int i3 = int(rndmPtr->flat() *acols.size());
      int acol3 = acols[i3];
      acols.erase(acols.begin() + i3);
      event.appendJunction(2,acol1,acol2,acol3);
    }

    // If we have more than two colour anti colour pairs
    // form junction anti junction pair.
    while (int(acols.size()) > 1) {
      int i1 = int(rndmPtr->flat() *cols.size());
      int col1 = cols[i1];
      cols.erase(cols.begin() + i1);
      int i2 = int(rndmPtr->flat() *cols.size());
      int col2 = cols[i2];
      cols.erase(cols.begin() + i2);
      int i3 = int(rndmPtr->flat() *acols.size());
      int acol1 = acols[i3];
      acols.erase(acols.begin() + i3);
      int i4 = int(rndmPtr->flat() *acols.size());
      int acol2 = acols[i4];
      acols.erase(acols.begin() + i4);
      int newCol = event.nextColTag();
      event.appendJunction(1, col1, col2, newCol);
      event.appendJunction(2, acol1, acol2, newCol);
    }

    // If we have one colour and one anti colour, form normal string.
    if (int(acols.size()) == 1) {
      int iCol = -1;
      for (int iPar = 0; iPar < event.size(); ++iPar)
        if (event[iPar].isFinal() && event[iPar].col() == cols[0])
          iCol = iPar;
      if (iCol == -1) {
        infoPtr->errorMsg("Warning in JunctionSplitting::SplitJunChain:"
          "Splitting multiple directly connected junctions failed");
        return false;
      }
      int iNew = event.copy(iCol, 76);
      event[iNew].col(acols[0]);
    }
  }

  // Delete the junctions from the event record.
  sort(junRem.begin(),junRem.end());
  reverse(junRem.begin(),junRem.end());
  for (int i = 0; i < int(junRem.size()); ++i)
    event.eraseJunction(junRem[i]);
  event.saveJunctionSize();

  return true;
}

//--------------------------------------------------------------------------

// Split junction pairs.
// If it has 3 connections everything fails.
// If it has 2 connections colapse into single string.
// If it has 1 connection, depend on the string length.

bool JunctionSplitting::splitJunPairs(Event& event,
  vector<vector< int > >& iPartonJun, vector<vector< int > >& iPartonAntiJun) {

  // Clear old memory.
  event.saveJunctionSize();
  vector<int> junRem;

  // Get junction chains.
  vector<vector<int> > junChains = colTrace.getJunChains(event);

  for (int i = 0; i < int(junChains.size()); ++i) {
    if (junChains[i].size() == 2) {
      vector<pair<int,int> > matchedLegs;
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          if (event.colJunction(junChains[i][0],j) ==
              event.colJunction(junChains[i][1],k))
            matchedLegs.push_back(make_pair(j,k));

      // Cannot handle 3 leg combination.
      if (matchedLegs.size() == 3) {
        infoPtr->errorMsg("Error in JunctionSplitting::SplitJunChain:"
          "Two junctions have all three legs connected.");
        return false;
      }

      // Split into string if two legs are combined.
      if (matchedLegs.size() == 2) {

        // Find first leg.
        int i1 = 0;
        if (matchedLegs[0].first != 1 && matchedLegs[1].first != 1) i1 = 1;
        if (matchedLegs[0].first != 2 && matchedLegs[1].first != 2) i1 = 2;

        // Find second leg.
        int j1 = 0;
        if (matchedLegs[0].second != 1 && matchedLegs[1].second != 1) j1 = 1;
        if (matchedLegs[0].second != 2 && matchedLegs[1].second != 2) j1 = 2;
        
        // Find corresponding colours.
        int col = event.colJunction(junChains[i][0],i1);
        int acol = event.colJunction(junChains[i][1],j1);
        if (event.kindJunction(junChains[i][1]) % 2 == 1)
          swap(col,acol);
        
        // Find index of anti particle.
        int iAcol = -1;
        for (int j = 0;j < event.size();++j)
          if (event[j].isFinal() && event[j].acol() == acol) {
            iAcol = j;
            break;
          }
        if (iAcol == -1) {
          infoPtr->errorMsg("Error in JunctionSplitting::SplitJunChain:"
            "Anti colour not found when combing two junctions to a string");
          return false;
        }
        
        // Update anti colour of anti particle.
        int iNew = event.copy(iAcol,66);
        event[iNew].acol(col);

        // Remove the junctions from the event record.
        junRem.push_back(junChains[i][0]);
        junRem.push_back(junChains[i][1]);
      }

      // Split into string if two legs are combined.
      if (matchedLegs.size() == 1) {

        // store junction numbers.
        int iJun = junChains[i][0];
        int iAnti = junChains[i][1];
        int iLeg = matchedLegs[0].first;
        int iAntiLeg = matchedLegs[0].second;
        if (event.kindJunction(iAnti) % 2 == 1)  {
          swap(iJun, iAnti);
          swap(iLeg, iAntiLeg);
        }

        // Find the junctions in the parton list.
        int iJunList = -1, iAntiList = -1;
        for (int l = 0;l < int(iPartonJun.size()); ++l)
          if (- iPartonJun[l][0]/10 - 1 == iJun) {
            iJunList = l;
            break;
          }
        
        for (int l = 0;l < int(iPartonAntiJun.size()); ++l)
          if (- iPartonAntiJun[l][0]/10 - 1 == iAnti) {
            iAntiList = l;
            break;
          }

        // Fill in vector of the legs content.
        vector<vector <int> > iJunLegs;
        iJunLegs.resize(3);
        int leg = -1;

        for (int l = 0; l < int(iPartonJun[iJunList].size()); ++l) {
          if ( iPartonJun[iJunList][l]/10 == iPartonJun[iJunList][0]/10) ++leg;
          iJunLegs[leg].push_back(iPartonJun[iJunList][l]);
        }

         // Fill in vector of the legs content.
        vector<vector <int> > iAntiLegs;
        iAntiLegs.resize(3);
        leg = -1;
        for (int l = 0; l < int(iPartonAntiJun[iAntiList].size()); ++l) {
          if ( iPartonAntiJun[iAntiList][l]/10
            == iPartonAntiJun[iAntiList][0]/10) ++leg;
          iAntiLegs[leg].push_back(iPartonAntiJun[iAntiList][l]);
        }

        // Identify the two external legs of either junction.
        vector<int>& iJunLeg0 = (iLeg == 0) ? iJunLegs[1] : iJunLegs[0];
        vector<int>& iJunLeg1 = (iLeg == 2) ? iJunLegs[1] : iJunLegs[2];
        vector<int>& iAntiLeg0 = (iAntiLeg == 0) ? iAntiLegs[1] : iAntiLegs[0];
        vector<int>& iAntiLeg1 = (iAntiLeg == 2) ? iAntiLegs[1] : iAntiLegs[2];

        // Simplified procedure: mainly study first parton on each leg.
        Vec4 pJunLeg0 = event[ iJunLeg0[1] ].p();
        Vec4 pJunLeg1 = event[ iJunLeg1[1] ].p();
        Vec4 pAntiLeg0 = event[ iAntiLeg0[1] ].p();
        Vec4 pAntiLeg1 = event[ iAntiLeg1[1] ].p();

      // Starting frame hopefully intermediate to two junction directions.
      Vec4 pStart = pJunLeg0 / pJunLeg0.e() + pJunLeg1 / pJunLeg1.e()
        + pAntiLeg0 / pAntiLeg0.e() + pAntiLeg1 / pAntiLeg1.e();

      // Loop over iteration to junction/antijunction rest frames (JRF/ARF).
      RotBstMatrix MtoJRF, MtoARF;
      Vec4 pInJRF[3], pInARF[3];
      for (int iJunLocal = 0; iJunLocal < 2; ++iJunLocal) {
        int offset = (iJunLocal == 0) ? 0 : 2;

        // Iterate from system rest frame towards the junction rest frame.
        RotBstMatrix MtoRF, Mstep;
        MtoRF.bstback(pStart);
        Vec4 pInRF[4];
        int iter = 0;
        do {
          ++iter;

          // Find rest-frame momenta on the three sides of the junction.
          // Only consider first parton on each leg, for simplicity.
          pInRF[0 + offset] = pJunLeg0;
          pInRF[1 + offset] = pJunLeg1;
          pInRF[2 - offset] = pAntiLeg0;
          pInRF[3 - offset] = pAntiLeg1;
          for (int l = 0; l < 4; ++l) pInRF[l].rotbst(MtoRF);

          // For third side add both legs beyond other junction, weighted.
          double wt2 = 1. - exp( -pInRF[2].e() / eNormJunction);
          double wt3 = 1. - exp( -pInRF[3].e() / eNormJunction);
          pInRF[2] = wt2 * pInRF[2] + wt3 * pInRF[3];

          // Find new junction rest frame from the set of momenta.
          Mstep = stringFrag.junctionRestFrame( pInRF[0], pInRF[1], pInRF[2]);
          MtoRF.rotbst( Mstep );
        } while (iter < 3 || (Mstep.deviation() > CONVJNREST
          && iter < NTRYJNREST) );

        // Store final boost and rest-frame (weighted) momenta.
        if (iJunLocal == 0) {
          MtoJRF = MtoRF;
          for (int l = 0; l < 3; ++l) pInJRF[l] = pInRF[l];
        } else {
          MtoARF = MtoRF;
          for (int l = 0; l < 3; ++l) pInARF[l] = pInRF[l];
        }
      }

      // Opposite operations: boost from JRF/ARF to original system.
      RotBstMatrix MfromJRF = MtoJRF;
      MfromJRF.invert();
      RotBstMatrix MfromARF = MtoARF;
      MfromARF.invert();

      // Velocity vectors of junctions and momentum of legs in lab frame.
      Vec4 vJun(0., 0., 0., 1.);
      vJun.rotbst(MfromJRF);
      Vec4 vAnti(0., 0., 0., 1.);
      vAnti.rotbst(MfromARF);
      Vec4 pLabJ[3], pLabA[3];
      for (int l = 0; l < 3; ++l) {
        pLabJ[l] = pInJRF[l];
        pLabJ[l].rotbst(MfromJRF);
        pLabA[l] = pInARF[l];
        pLabA[l].rotbst(MfromARF);
      }

      // Calculate Lambda-measure length of three possible topologies.
      double vJvA = vJun * vAnti;
      double vJvAe2y = vJvA + sqrt(vJvA*vJvA - 1.);
      double lambdaJA = stringLength.getJuncLength(pInJRF[0], pInJRF[1],
        pInARF[0], pInARF[1]);

      double lambda00 = stringLength.getStringLength(pLabJ[0], pLabA[0]) +
        stringLength.getStringLength(pLabJ[1], pLabA[1]);

      double lambda01 = stringLength.getStringLength(pLabJ[0], pLabA[1]) +
        stringLength.getStringLength(pLabJ[1], pLabA[0]);

      // Case when either topology without junctions is the shorter one.
      if (lambdaJA > min( lambda00, lambda01) && allowDoubleJunRem) {

        // Find indices of particles.
        int iCol1 = iJunLeg0[1];
        int iCol2 = iJunLeg1[1];
        int iAcol1 = iAntiLeg0[1];
        int iAcol2 = iAntiLeg1[1];
        if (lambda00 > lambda01)
          swap(iAcol1, iAcol2);

        // Change the colour index and mark junctions to be removed.
        int iNew1 = event.copy(iAcol1, 66);
        event[iNew1].acol(event[iCol1].col());
        int iNew2 = event.copy(iAcol2, 66);
        event[iNew2].acol(event[iCol2].col());
        junRem.push_back(junChains[i][0]);
        junRem.push_back(junChains[i][1]);
        continue;
      }

      // Case where junction and antijunction to be separated.
      // Shuffle (p+/p-)  momentum of order <mThad> between systems,
      // times 2/3 for 120 degree in JRF, times 1/2 for two legs,
      // but not more than half of what nearest parton carries.
      double eShift = MTHAD / (3. * sqrt(vJvAe2y));
      double fracJ0 = min(0.5, eShift / pInJRF[0].e());
      double fracJ1 = min(0.5, eShift / pInJRF[1].e());
      Vec4 pFromJun = fracJ0 * pJunLeg0 + fracJ1 * pJunLeg1;
      double fracA0 = min(0.5, eShift / pInARF[0].e());
      double fracA1 = min(0.5, eShift / pInARF[1].e());
      Vec4 pFromAnti = fracA0 * pAntiLeg0 + fracA1 * pAntiLeg1;

      // Pick a new quark at random; for simplicity no diquarks.
      int idQ = flavSel.pickLightQ();

      // Copy junction partons with scaled-down momenta and update legs.
      int mother1 = min(iJunLeg0[1], iJunLeg1[1]);
      int mother2 = max(iJunLeg0[1], iJunLeg1[1]);
      int iNew1 = event.copy(iJunLeg0[1], 76);
      event[iNew1].rescale5(1. - fracJ0);
      iJunLeg0[1] = iNew1;
      int iNew2 = event.copy(iJunLeg1[1], 76);
      event[iNew2].rescale5(1. - fracJ1);
      iJunLeg1[1] = iNew2;

      // Update junction colour and store quark with antijunction momentum.
      // Store history as 2 -> 3  step for consistency.
      int colQ = event.nextColTag();
      event.endColJunction(iJun, iLeg, colQ);
      event.colJunction(iJun, iLeg, colQ);
      int iNewJ = event.append( idQ, 76, mother1, mother2, 0, 0,
        colQ, 0, pFromAnti, pFromAnti.mCalc() );
      event[mother1].daughters( iNew1, iNewJ);
      event[mother2].daughters( iNew1, iNewJ);
      event[iNew1].mothers( mother1, mother2);
      event[iNew2].mothers( mother1, mother2);

      // Copy anti junction partons with scaled-down momenta and update legs.
      mother1 = min(iAntiLeg0[1], iAntiLeg1[1]);
      mother2 = max(iAntiLeg0[1], iAntiLeg1[1]);
      iNew1 = event.copy(iAntiLeg0[1], 76);
      event[iNew1].rescale5(1. - fracA0);
      iAntiLeg0[1] = iNew1;
      iNew2 = event.copy(iAntiLeg1[1], 76);
      event[iNew2].rescale5(1. - fracA1);
      iAntiLeg1[1] = iNew2;

      // Update antijunction anticolour and store antiquark with junction
      // momentum. Store history as 2 -> 3  step for consistency.
      int acolQ = event.nextColTag();
      event.endColJunction(iAnti, iAntiLeg, acolQ);
      event.colJunction(iAnti, iAntiLeg, acolQ);
      int iNewA = event.append( -idQ, 76, mother1, mother2, 0, 0,
        0, acolQ, pFromJun, pFromJun.mCalc() );
      event[mother1].daughters( iNew1, iNewA);
      event[mother2].daughters( iNew1, iNewA);
      event[iNew1].mothers( mother1, mother2);
      event[iNew2].mothers( mother1, mother2);

      // Done with splitting junction from antijunction.
      }
    }
  }

  // Delete the junctions from the event record.
  sort(junRem.begin(),junRem.end());
  reverse(junRem.begin(),junRem.end());
  for (int i = 0; i < int(junRem.size()); ++i)
    event.eraseJunction(junRem[i]);
  event.saveJunctionSize();

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Get the list of partons connected to the junctions.

bool JunctionSplitting::getPartonLists(Event& event,
  vector<vector< int > > & iPartonJun, vector<vector<int > >& iPartonAntiJun) {

  // Need to try and split junction structures.
  colTrace.setupColList(event);
  vector<int> iParton;
  iPartonJun.clear();
  iPartonAntiJun.clear();

  // Loop over junctions and collect all junctions.
  // Then afterwards collect all anti junctions.
  // This ensures that all gluons are collected on the junctions.
  for (int iJun = 0; iJun < event.sizeJunction(); ++iJun)
  if (event.remainsJunction(iJun)) {

    int kindJun = event.kindJunction(iJun);
    iParton.resize(0);

    if (kindJun % 2 != 1) continue;

    // Loop over junction legs.
    for (int iCol = 0; iCol < 3; ++iCol) {
      int indxCol = event.colJunction(iJun, iCol);
      iParton.push_back( -(10 + 10 * iJun + iCol) );
      // Junctions: find color ends.
      if (kindJun % 2 == 1 && !colTrace.traceFromAcol(indxCol,
       event, iJun, iCol, iParton))
        return false;
    }

    // Store the anti junction and junction list.
    int nNeg = 0;
    for (int i = 0; i < int(iParton.size()); ++i) if (iParton[i] < 0)
      ++nNeg;
    if (nNeg > 3 )
      iPartonJun.push_back(iParton);
  }

  // Loop over all anti junctions.
  for (int iJun = 0; iJun < event.sizeJunction(); ++iJun)
  if (event.remainsJunction(iJun)) {

    int kindJun = event.kindJunction(iJun);
    iParton.resize(0);

    if (kindJun % 2 != 0)
      continue;

    // Loop over junction legs.
    for (int iCol = 0; iCol < 3; ++iCol) {
      int indxCol = event.colJunction(iJun, iCol);
      iParton.push_back( -(10 + 10 * iJun + iCol) );
      // Antijunctions: find anticolor ends.
      if (kindJun % 2 == 0 && !colTrace.traceFromCol(indxCol,
       event, iJun, iCol, iParton))
        return false;
    }

    // Store the anti junction and junction list.
    int nNeg = 0;
    for (int i = 0; i < int(iParton.size()); ++i) if (iParton[i] < 0)
      ++nNeg;
    if (nNeg > 3 )
      iPartonAntiJun.push_back(iParton);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Change the anticolour of the particle that has acol to be col.

bool JunctionSplitting::setAcol(Event& event, int col, int acol) {

  // Update anticolour if it belongs to a particle.
  for (int j = 0;j < event.size(); ++j)
    if (event[j].isFinal() && event[j].acol() == acol) {
      int iNew = event.copy(j,66);
      event[iNew].acol(col);
      return true;
    }
  // Check if anti junction is connected to a junction.
  for (int j = 0;j < event.sizeJunction(); ++j)
    for (int jLeg = 0;jLeg < 3; ++jLeg)
      if (event.colJunction(j, jLeg) == acol) {
        event.colJunction(j, jLeg, col);
        return true;
      }

  // If no acol was found something went wrong.
  infoPtr->errorMsg("Error in JunctionSplitting::setAcol:"
     "Anti colour not found when combing two junctions to a string");
  return false;
}

//==========================================================================

} // end namespace Pythia8


