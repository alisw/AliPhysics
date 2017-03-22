// RHadrons.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the RHadrons class.

#include "RHadrons.h"

namespace Pythia8 {
 
//==========================================================================

// The RHadrons class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

const int RHadrons::IDRHADSB[14] = {  1000512, 1000522, 1000532,
  1000542, 1000552, 1005113, 1005211, 1005213, 1005223, 1005311,
  1005313, 1005321, 1005323, 1005333 };

const int RHadrons::IDRHADST[14] = {  1000612, 1000622, 1000632,
  1000642, 1000652, 1006113, 1006211, 1006213, 1006223, 1006311,
  1006313, 1006321, 1006323, 1006333 };

// Gluino code and list of gluino R-hadron codes.
const int RHadrons::IDRHADGO[38] = {  1000993, 1009113, 1009213, 
  1009223, 1009313, 1009323, 1009333, 1009413, 1009423, 1009433, 
  1009443, 1009513, 1009523, 1009533, 1009543, 1009553, 1091114, 
  1092114, 1092214, 1092224, 1093114, 1093214, 1093224, 1093314, 
  1093324, 1093334, 1094114, 1094214, 1094224, 1094314, 1094324, 
  1094334, 1095114, 1095214, 1095224, 1095314, 1095324, 1095334 };

// Allow a few tries for flavour selection since failure is allowed.
const int RHadrons::NTRYMAX = 10;

// Safety margin (in GeV) when constructing kinematics of system.
const double RHadrons::MSAFETY = 0.1;

// Maximal energy to borrow for gluon to insert on leg in to junction.
const double RHadrons::EGBORROWMAX = 4.;

//--------------------------------------------------------------------------

// Main routine to initialize R-hadron handling.

bool RHadrons::init( Info* infoPtrIn, Settings& settings, 
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn) {

  // Store input pointers for future use. 
  infoPtr          = infoPtrIn;
  particleDataPtr  = particleDataPtrIn;  
  rndmPtr          = rndmPtrIn;

  // Flags and parameters related to R-hadron formation and decay.
  allowRH          = settings.flag("RHadrons:allow");
  maxWidthRH       = settings.parm("RHadrons:maxWidth");
  idRSb            = settings.mode("RHadrons:idSbottom");
  idRSt            = settings.mode("RHadrons:idStop");
  idRGo            = settings.mode("RHadrons:idGluino");
  setMassesRH      = settings.flag("RHadrons:setMasses");
  probGluinoballRH = settings.parm("RHadrons:probGluinoball");
  mOffsetCloudRH   = settings.parm("RHadrons:mOffsetCloud");
  mCollapseRH      = settings.parm("RHadrons:mCollapse");
  diquarkSpin1RH   = settings.parm("RHadrons:diquarkSpin1"); 

  // Check whether sbottom, stop or gluino should form R-hadrons. 
  allowRSb         = allowRH && idRSb > 0 
    && (particleDataPtr->mWidth(idRSb) < maxWidthRH);
  allowRSt         = allowRH && idRSt > 0 
    && (particleDataPtr->mWidth(idRSt) < maxWidthRH);
  allowRGo         = allowRH && idRGo > 0 
    && (particleDataPtr->mWidth(idRGo) < maxWidthRH);
  allowSomeR       = allowRSb || allowRSt || allowRGo;

  // Set nominal masses of sbottom R-mesons and R-baryons.
  if (allowRSb) {
    m0Sb = particleDataPtr->m0(idRSb);
    if (setMassesRH) {
      for (int i = 0; i < 14; ++i) {
        int idR = IDRHADSB[i]; 
        double m0RHad = m0Sb + mOffsetCloudRH;
        m0RHad += particleDataPtr->constituentMass( (idR%100)/10);
        if (i > 4) 
        m0RHad += particleDataPtr->constituentMass( (idR%1000)/100 );
        particleDataPtr->m0( idR, m0RHad);    
      }
    }

    // Set widths and lifetimes of sbottom R-hadrons.
    double mWidthRHad = particleDataPtr->mWidth(idRSb);
    double tau0RHad   = particleDataPtr->tau0(  idRSb);
    for (int i = 0; i < 14; ++i) {
      particleDataPtr->mWidth( IDRHADSB[i], mWidthRHad);
      particleDataPtr->tau0(   IDRHADSB[i],   tau0RHad);
    }
  }

  // Set nominal masses of stop R-mesons and R-baryons.
  if (allowRSt) {
    m0St = particleDataPtr->m0(idRSt);
    if (setMassesRH) {
      for (int i = 0; i < 14; ++i) {
        int idR = IDRHADST[i]; 
        double m0RHad = m0St + mOffsetCloudRH;
        m0RHad += particleDataPtr->constituentMass( (idR%100)/10);
        if (i > 4) 
        m0RHad += particleDataPtr->constituentMass( (idR%1000)/100 );
        particleDataPtr->m0( idR, m0RHad);    
      }
    }

    // Set widths and lifetimes of stop R-hadrons.
    double mWidthRHad = particleDataPtr->mWidth(idRSt);
    double tau0RHad   = particleDataPtr->tau0(  idRSt);
    for (int i = 0; i < 14; ++i) {
      particleDataPtr->mWidth( IDRHADST[i], mWidthRHad);
      particleDataPtr->tau0(   IDRHADST[i],   tau0RHad);
    }
  }

  // Set nominal masses of gluino R-glueballs, R-mesons and R-baryons.
  if (allowRGo) {
    m0Go = particleDataPtr->m0(idRGo);
    if (setMassesRH) {
      particleDataPtr->m0( IDRHADGO[0], m0Go + 2. * mOffsetCloudRH 
        + particleDataPtr->constituentMass(21) );
      for (int i = 1; i < 38; ++i) {
        int idR = IDRHADGO[i]; 
        double m0RHad = m0Go + 2. * mOffsetCloudRH;
        m0RHad += particleDataPtr->constituentMass( (idR%1000)/100 );
        m0RHad += particleDataPtr->constituentMass( (idR%100)/10);
        if (i > 15) 
        m0RHad += particleDataPtr->constituentMass( (idR%10000)/1000 );
        particleDataPtr->m0( idR, m0RHad);    
      }
    }

    // Set widths and lifetimes of gluino R-hadrons.
    double mWidthRHad = particleDataPtr->mWidth(idRGo);
    double tau0RHad   = particleDataPtr->tau0(  idRGo);
    for (int i = 0; i < 38; ++i) {
      particleDataPtr->mWidth( IDRHADGO[i], mWidthRHad);
      particleDataPtr->tau0(   IDRHADGO[i],   tau0RHad);
    }
  }   

  // Done. 
  return true;

}
//--------------------------------------------------------------------------

// Check if a given particle can produce R-hadrons.

bool RHadrons::givesRHadron( int id) {
  if (allowRSb && abs(id) == idRSb) return true;
  if (allowRSt && abs(id) == idRSt) return true;
  if (allowRGo && id == idRGo) return true;
  return false;

}

//--------------------------------------------------------------------------

// Produce R-hadrons by fragmenting them off from existing strings.

bool RHadrons::produce( ColConfig& colConfig, Event& event) {

  // Check whether some sparticles are allowed to hadronize.
  if (!allowSomeR) return true;

  // Reset arrays and values.
  iBefRHad.resize(0);
  iCreRHad.resize(0);
  iRHadron.resize(0);
  iAftRHad.resize(0);
  isTriplet.resize(0);
  nRHad = 0;

  // Loop over event and identify hadronizing sparticles.
  for (int i = 0; i < event.size(); ++i) 
   if (event[i].isFinal() && givesRHadron(event[i].id())) { 
    iBefRHad.push_back(i);
    iCreRHad.push_back(i);
    iRHadron.push_back(0);
    iAftRHad.push_back(0);
    isTriplet.push_back(true);
  } 
  nRHad = iRHadron.size();
  
  // Done if no hadronizing sparticles.
  if (nRHad == 0) return true;

  // Max two R-hadrons. Randomize order of processing.
  if (nRHad > 2) {
     infoPtr->errorMsg("Error in RHadrons::produce: "
       "cannot handle more than two R-hadrons");
     return false;
  }
  if (nRHad > 1 && rndmPtr->flat() > 0.5) swap( iBefRHad[0], iBefRHad[1]);

  // Split a system with both a sparticle and a junction.
  iBef = iBefRHad[0];  
  iSys = colConfig.findSinglet( iBef);
  systemPtr = &colConfig[iSys];
  if (systemPtr->hasJunction && !splitOffJunction( colConfig, event)) {
    infoPtr->errorMsg("Error in RHadrons::produce: "
      "cannot handle system with junction");
    return false;
  }
  if (nRHad == 2) {
    iBef = iBefRHad[1];  
    iSys = colConfig.findSinglet( iBefRHad[1]);
    systemPtr = &colConfig[iSys];
    if (systemPtr->hasJunction && !splitOffJunction( colConfig, event)) {
      infoPtr->errorMsg("Error in RHadrons::produce: "
        "cannot handle system with junction");
      return false;
    }
  }

  // Open up a closed gluon/gluino loop.
  iBef = iBefRHad[0];  
  iSys = colConfig.findSinglet( iBef);
  systemPtr = &colConfig[iSys];
  if (systemPtr->isClosed && !openClosedLoop( colConfig, event)) {
    infoPtr->errorMsg("Error in RHadrons::produce: "
      "cannot open up closed gluon/gluino loop");
    return false;
  }
  if (nRHad == 2) {
    iBef = iBefRHad[1];  
    iSys = colConfig.findSinglet( iBefRHad[1]);
    systemPtr = &colConfig[iSys];
    if (systemPtr->isClosed && !openClosedLoop( colConfig, event)) {
      infoPtr->errorMsg("Error in RHadrons::produce: "
        "cannot open up closed gluon/gluino loop");
      return false;
    }
  }

  // Split up a colour singlet system that should give two R-hadrons.
  if (nRHad == 2) {
    int iSys1 = colConfig.findSinglet( iBefRHad[0]);
    int iSys2 = colConfig.findSinglet( iBefRHad[1]);
    if (iSys2 == iSys1) { 
      iSys = iSys1;
      systemPtr = &colConfig[iSys];
      if ( !splitSystem( colConfig, event) ) { 
        infoPtr->errorMsg("Error in RHadrons::produce: "
          "failed to handle two sparticles in same system");
        return false;
      }
    } 
  }
    
  // Loop over R-hadrons to be formed. Find its colour singlet system.
  for (iRHad = 0; iRHad < nRHad; ++iRHad) {
    iBef = iBefRHad[iRHad];  
    iSys = colConfig.findSinglet( iBef);
    if (iSys < 0) {
      infoPtr->errorMsg("Error in RHadrons::produce: "
        "sparticle not in any colour singlet");
      return false;
    }
    systemPtr = &colConfig[iSys];

    // For now don't handle systems involving junctions or loops.
    if (systemPtr->hasJunction) {
      infoPtr->errorMsg("Error in RHadrons::produce: "
        "cannot handle system with junction");
      return false;
    }
    if (systemPtr->isClosed) {
      infoPtr->errorMsg("Error in RHadrons::produce: "
        "cannot handle closed colour loop");
      return false;
    }

    // Handle formation of R-hadron separately for gluino and squark.
    if (event[iBef].id() == idRGo) isTriplet[iRHad] = false;
    bool formed = (isTriplet[iRHad]) ? produceSquark( colConfig, event)
                                     : produceGluino( colConfig, event);
    if (!formed) return false;

  // End of loop over R-hadrons. Done.
  } 
  return true;

}

//--------------------------------------------------------------------------

// Decay R-hadrons by resolving them into string systems and letting the
// heavy unstable particle decay as normal.

bool RHadrons::decay( Event& event) {

  // Loop over R-hadrons to decay. 
  for (iRHad = 0; iRHad < nRHad; ++iRHad) {
    int    iRNow  = iRHadron[iRHad]; 
    int    iRBef  = iBefRHad[iRHad];
    int    idRHad = event[iRNow].id();
    double mRHad  = event[iRNow].m();
    double mRBef  = event[iRBef].m();
    int    iR0    = 0;
    int    iR2    = 0; 

    // Find flavour content of squark or gluino R-hadron.
    pair<int,int> idPair = (isTriplet[iRHad]) 
      ? fromIdWithSquark( idRHad) : fromIdWithGluino( idRHad);
    int id1 = idPair.first;
    int id2 = idPair.second;

    // Sharing of momentum: the squark/gluino should be restored
    // to original mass, but error if negative-mass spectators.
    double fracR = mRBef / mRHad;
    if (fracR >= 1.) {
      infoPtr->errorMsg("Error in RHadrons::decay: "
          "too low R-hadron mass for decay");
      return false;
    }

    // Squark: new colour needed in the breakup.
    if (isTriplet[iRHad]) {
      int colNew = event.nextColTag();
      int col    = (event[iRBef].col() != 0) ? colNew : 0;
      int acol   = (col == 0) ? colNew : 0; 
      
      // Store the constituents of a squark R-hadron.
      iR0 = event.append( id1, 106, iRNow, 0, 0, 0, col, acol,
        fracR * event[iRNow].p(), fracR * mRHad, 0.);
      iR2 = event.append( id2, 106, iRNow, 0, 0, 0, acol, col, 
        (1. - fracR) * event[iRNow].p(), (1. - fracR) * mRHad, 0.);

    // Gluino: set mass sharing between two spectators.
    } else {
      double m1Eff  = particleDataPtr->constituentMass(id1) + mOffsetCloudRH;
      double m2Eff  = particleDataPtr->constituentMass(id2) + mOffsetCloudRH;
      double frac1 = (1. - fracR) * m1Eff / ( m1Eff + m2Eff); 
      double frac2 = (1. - fracR) * m2Eff / ( m1Eff + m2Eff); 
   
      // Two new colours needed in the breakups.
      int col1 = event.nextColTag();
      int col2 = event.nextColTag();

      // Store the constituents of a gluino R-hadron.
      iR0 = event.append( idRGo, 106, iRNow, 0, 0, 0, col2, col1,
        fracR * event[iRNow].p(), fracR * mRHad, 0.);
      event.append( id1, 106, iRNow, 0, 0, 0, col1, 0, 
        frac1 * event[iRNow].p(), frac1 * mRHad, 0.);
      iR2 = event.append( id2, 106, iRNow, 0, 0, 0, 0, col2, 
        frac2 * event[iRNow].p(), frac2 * mRHad, 0.);
    }

    // Mark R-hadron as decayed and update history.
    event[iRNow].statusNeg();
    event[iRNow].daughters( iR0, iR2);
    iAftRHad[iRHad] = iR0;

    // Set secondary vertex for decay products, but no lifetime.
    Vec4 vDec = event[iRNow].vProd() + event[iRNow].tau()
              * event[iR0].p() / event[iR0].m();
    for (int iRd = iR0; iRd <= iR2; ++iRd) event[iRd].vProd( vDec);

  // End loop over R-hadron decays, based on velocity of squark.
  
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Split a system that contains both a sparticle and a junction. 

bool RHadrons::splitOffJunction( ColConfig& colConfig, Event& event) {

  // Classify system into three legs, and find sparticle location.
  vector<int> leg1, leg2, leg3;
  int iLegSP = 0;
  int iPosSP = 0;
  int iLeg = 0;
  int iPos = 0;
  for (int i = 0; i < systemPtr->size(); ++i) {
    ++iPos;
    int iP = systemPtr->iParton[i];
    if (iP < 0) {
      ++iLeg;
      iPos = -1;
    } else if (iLeg == 1) leg1.push_back( iP);
    else if   (iLeg == 2) leg2.push_back( iP);
    else if   (iLeg == 3) leg3.push_back( iP);
    if (iP == iBef) {
      iLegSP = iLeg;
      iPosSP = iPos;
    }
  }
  if (iLegSP == 0) return false;
  
  // Swap so leg 1 contains sparticle. If not innermost sparticle then
  // skip for now, and wait for this other (gluino!) to be split off. 
  if      (iLegSP == 2) swap(leg2, leg1);
  else if (iLegSP == 3) swap(leg3, leg1); 
  for (int i = 0; i < iPosSP; ++i)
    if (event[leg1[i]].id() != 21) return true;
 
  // No gluon between sparticle and junction: find kinetic energy of system.
  if (iPosSP == 0) {
    Vec4 pSP  = event[iBef].p();
    Vec4 pRec = 0.;
    for (int i = 0; i < int(leg2.size()); ++i) pRec += event[leg2[i]].p();
    for (int i = 0; i < int(leg3.size()); ++i) pRec += event[leg3[i]].p();
    double mSys  = (pSP + pRec).mCalc();
    double mSP   = pSP.mCalc();
    double mRec  = pRec.mCalc();
    double eKin  = mSys - mSP - mRec;

    // Insert new gluon that borrows part of kinetic energy.
    double mNewG  = EGBORROWMAX * eKin / (EGBORROWMAX + eKin) ;
    Vec4   pNewG  = (mNewG / mSys) * (pSP + pRec);
    int    newCol = event.nextColTag();
    bool   isCol  = (event[leg1.back()].col() > 0);
    int    colNG  = (isCol)? newCol :  event[iBef].acol();
    int    acolNG = (isCol) ? event[iBef].col() : newCol;
    int    iNewG  = event.append( 21, 101, iBef, 0, 0, 0, colNG, acolNG, 
      pNewG, mNewG, 0.);
    leg1.insert( leg1.begin(), iNewG);
    ++iPosSP;

    // Boosts for remainder systems that gave up energy.
    double mNewSys = mSys - mNewG;
    double pAbsOld = 0.5 * sqrtpos( pow2(mSys*mSys - mSP*mSP - mRec*mRec)
                   - pow2(2. * mSP * mRec) ) / mSys;
    double pAbsNew = 0.5 * sqrtpos( pow2(mNewSys*mNewSys - mSP*mSP - mRec*mRec)
                   - pow2(2. * mSP * mRec) ) / mNewSys;
    RotBstMatrix MtoCM, MfromCM, MSP, MRec;
    MtoCM.toCMframe( pSP, pRec);
    MfromCM = MtoCM;
    MfromCM.invert(); 
    MSP = MtoCM;
    MSP.bst( 0., 0., -pAbsOld / sqrt(mSP * mSP + pAbsOld * pAbsOld));
    MSP.bst( 0., 0.,  pAbsNew / sqrt(mSP * mSP + pAbsNew * pAbsNew)); 
    MSP.rotbst( MfromCM);
    MRec = MtoCM;
    MRec.bst( 0., 0.,  pAbsOld / sqrt(mRec * mRec + pAbsOld * pAbsOld));
    MRec.bst( 0., 0., -pAbsNew / sqrt(mRec * mRec + pAbsNew * pAbsNew)); 
    MRec.rotbst( MfromCM);

    // Copy down recoling partons and boost their momenta.
    int iNewSP  = event.copy( iBef, 101);
    event[iNewSP].mother2(0);
    event[iBef].daughter1(iNewG);
    event[iNewSP].rotbst( MSP);
    leg1[iPosSP]   = iNewSP;
    if (iBefRHad[0] == iBef) iBefRHad[0] = iNewSP;
    else if (nRHad > 1 && iBefRHad[1] == iBef) iBefRHad[1] = iNewSP;
    iBef = iNewSP;
    for (int i = 0; i < int(leg2.size()); ++i) {
      int iCopy = event.copy( leg2[i], 101);  
      event[iCopy].rotbst( MRec);
      if (iBefRHad[0] == leg2[i]) iBefRHad[0] = iCopy;
      else if (nRHad > 1 && iBefRHad[1] == leg2[i]) iBefRHad[1] = iCopy;
      leg2[i] = iCopy;
    }   
    for (int i = 0; i < int(leg3.size()); ++i) {
      int iCopy = event.copy( leg3[i], 101);  
      event[iCopy].rotbst( MRec);
      if (iBefRHad[0] == leg3[i]) iBefRHad[0] = iCopy;
      else if (nRHad > 1 && iBefRHad[1] == leg3[i]) iBefRHad[1] = iCopy;
      leg3[i]   = iCopy;
    }
 
  // Now always at least one gluon between sparticle and junction.  
  }

  // Find gluon with largest energy in sparticle rest frame.
  int    iGspl = 0;
  double eGspl = event[leg1[0]].p() * event[iBef].p();
  for (int i = 1; i < iPosSP; ++i) {
    double eGnow = event[leg1[i]].p() * event[iBef].p();
    if (eGnow > eGspl) {
      iGspl = i;
      eGspl = eGnow;
    }
  }
  int iG = leg1[iGspl];
   
  // Split this gluon into a collinear quark.antiquark pair.
  int idNewQ = flavSelPtr->pickLightQ(); 
  int iNewQ  = event.append(  idNewQ, 101, iG, 0, 0, 0, event[iG].col(), 0, 
    0.5 * event[iG].p(), 0.5 * event[iG].m(), 0.);
  int iNewQb = event.append( -idNewQ, 101, iG, 0, 0, 0, 0, event[iG].acol(), 
    0.5 * event[iG].p(), 0.5 * event[iG].m(), 0.);
  event[iG].statusNeg();
  event[iG].daughters( iNewQ, iNewQb); 
  if (event[leg1.back()].col() == 0) swap( iNewQ, iNewQb);

  // Collect two new systems after split.
  vector<int> iNewSys1, iNewSys2;
  iNewSys1.push_back( iNewQb);
  for (int i = iGspl + 1; i < int(leg1.size()); ++i)
    iNewSys1.push_back( leg1[i]);
  iNewSys2.push_back( -10);
  for (int i = 0; i < iGspl; ++i) iNewSys2.push_back( leg1[i]);
  iNewSys2.push_back( iNewQ);
  iNewSys2.push_back( -11);
  for (int i = 0; i < int(leg2.size()); ++i) iNewSys2.push_back( leg2[i]);
  iNewSys2.push_back( -12);
  for (int i = 0; i < int(leg3.size()); ++i) iNewSys2.push_back( leg3[i]);

  // Remove old system and insert two new ones.
  colConfig.erase(iSys);
  colConfig.insert( iNewSys1, event);
  colConfig.insert( iNewSys2, event);   

  // Done. 
  return true;

}

//--------------------------------------------------------------------------

// Open up a closed gluon/gluino loop.
  
bool RHadrons::openClosedLoop( ColConfig& colConfig, Event& event) {

  // Find gluon with largest energy in gluino rest frame.
  int    iGspl = -1;
  double eGspl = 0.;
  for (int i = 0; i < systemPtr->size(); ++i) {
    int  iGNow = systemPtr->iParton[i];
    if (event[iGNow].id() == 21) {
      double eGnow = event[iGNow].p() * event[iBef].p();
      if (eGnow > eGspl) {
        iGspl = i;
        eGspl = eGnow;
      }
    }
  }
  if (iGspl == -1) return false;
   
  // Split this gluon into a collinear quark.antiquark pair.
  int iG     = systemPtr->iParton[iGspl];
  int idNewQ = flavSelPtr->pickLightQ(); 
  int iNewQ  = event.append(  idNewQ, 101, iG, 0, 0, 0, event[iG].col(), 0, 
    0.5 * event[iG].p(), 0.5 * event[iG].m(), 0.);
  int iNewQb = event.append( -idNewQ, 101, iG, 0, 0, 0, 0, event[iG].acol(), 
    0.5 * event[iG].p(), 0.5 * event[iG].m(), 0.);
  event[iG].statusNeg();
  event[iG].daughters( iNewQ, iNewQb); 
   
  // Order partons in new system. Note order of colour flow.
  int iNext = iGspl + 1;
  if (iNext == systemPtr->size()) iNext = 0; 
  if (event[ systemPtr->iParton[iNext]].acol() != event[iNewQ].col())
    swap( iNewQ, iNewQb);        
  vector<int> iNewSys;
  iNewSys.push_back( iNewQ);
  for (int i = iGspl + 1; i < systemPtr->size(); ++i)
    iNewSys.push_back( systemPtr->iParton[i]);
  for (int i = 0; i < iGspl; ++i)
    iNewSys.push_back( systemPtr->iParton[i]);
  iNewSys.push_back( iNewQb);  

  // Erase the old system and insert the new one instead.
  colConfig.erase(iSys);
  colConfig.insert( iNewSys, event);

  // Done. 
  return true;

}

//--------------------------------------------------------------------------

// Split a single colour singlet that contains two sparticles.
// To fix : if nLeg >= 3 && mMin large handle as in nLeg == 1??

bool RHadrons::splitSystem( ColConfig& colConfig, Event& event) {

  // First and second R-hadron mother. Number of legs in between.
  int iFirst  = -1;
  int iSecond = -1;
  for (int i = 0; i < int(systemPtr->size()); ++i) {
    int  iTmp   = systemPtr->iParton[i];
    if ( givesRHadron( event[iTmp].id()) ) { 
      if (iFirst == -1) iFirst  = i;
      else              iSecond = i;
    }
  }
  int nLeg = iSecond - iFirst;

  // New flavour pair for breaking the string, and its mass.
  int    idNewQ = flavSelPtr->pickLightQ();
  double mNewQ  = particleDataPtr->constituentMass( idNewQ);
  vector<int> iNewSys1, iNewSys2;

  // If sparticles are neighbours then need new q-qbar pair in between.
  if (nLeg == 1) {

    // The location of the two sparticles.
    int i1Old = systemPtr->iParton[iFirst];
    int i2Old = systemPtr->iParton[iSecond];

    // Take a fraction of momentum to give breakup quark pair.
    double mHat = (event[i1Old].p() + event[i2Old].p()).mCalc();
    double mMax = mHat - event[i1Old].m() - event[i2Old].m(); 
    if (mMax < 2. * (mNewQ + MSAFETY)) return false;
    double mEff = min( 2. * (mNewQ + mOffsetCloudRH), mMax - 2. * MSAFETY);
    double frac = mEff / mHat;
    Vec4   pEff = frac * (event[i1Old].p() + event[i2Old].p());
  
    // New kinematics by (1) go to same mHat with bigger masses, and 
    // (2) rescale back down to original masses and reduced mHat.
    Vec4 p1New, p2New;
    if ( !newKin( event[i1Old].p(), event[i2Old].p(), 
      event[i1Old].m() / (1. - frac), event[i2Old].m() / (1. - frac), 
      p1New, p2New) ) return false; 
    p1New *= 1. - frac;
    p2New *= 1. - frac;

    // Fill in new partons after branching, with correct colour/flavour sign.
    int i1New, i2New, i3New, i4New;
    int newCol = event.nextColTag();
    i1New = event.copy( i1Old, 101);
    if (event[i2Old].acol() == event[i1Old].col()) {
      i3New = event.append( -idNewQ, 101, i1Old, 0, 0, 0, 
        0, event[i2Old].acol(), 0.5 * pEff, 0.5 * mEff, 0.);
      i2New = event.copy( i2Old, 101);
      event[i2New].acol( newCol);
      i4New = event.append(  idNewQ, 101, i2Old, 0, 0, 0, 
        newCol, 0, 0.5 * pEff, 0.5 * mEff, 0.); 
    } else {
      i3New = event.append(  idNewQ, 101, i1Old, 0, 0, 0, 
        event[i2Old].col(), 0, 0.5 * pEff, 0.5 * mEff, 0.);
      i2New = event.copy( i2Old, 101);
      event[i2New].col( newCol);
      i4New = event.append( -idNewQ, 101, i2Old, 0, 0, 0, 
        0, newCol, 0.5 * pEff, 0.5 * mEff, 0.); 
    }

    // Modify momenta and history. For iBotCopyId tracing asymmetric 
    // bookkeeping where one sparticle is radiator and the other recoiler.
    event[i1New].p( p1New);
    event[i2New].p( p2New);
    event[i1Old].daughters( i1New, i3New);
    event[i1New].mother2( 0);
    event[i2Old].daughters( i2New, i4New);
    event[i2New].mother2( 0);
    iBefRHad[0] = i1New;
    iBefRHad[1] = i2New;
 
    // Partons in the two new systems.
    for (int i = 0; i < iFirst; ++i) 
      iNewSys1.push_back( systemPtr->iParton[i]);
    iNewSys1.push_back( i1New);
    iNewSys1.push_back( i3New);
    iNewSys2.push_back( i4New);
    iNewSys2.push_back( i2New);
    for (int i = iSecond + 1; i < int(systemPtr->size()); ++i) 
      iNewSys2.push_back( systemPtr->iParton[i]);

  // If one gluon between sparticles then split it into a
  // collinear q-qbar pair.
  } else if (nLeg == 2) {

    // Gluon to be split and its two daughters.
    int iOld  = systemPtr->iParton[iFirst + 1];
    int i1New = event.append(  idNewQ, 101, iOld, 0, 0, 0, 
      event[iOld].col(), 0, 0.5 * event[iOld].p(), 
      0.5 * event[iOld].m(), 0.);
    int i2New = event.append( -idNewQ, 101, iOld, 0, 0, 0, 
      0, event[iOld].acol(), 0.5 * event[iOld].p(), 
      0.5 * event[iOld].m(), 0.);
    event[iOld].statusNeg();
    event[iOld].daughters( i1New, i2New);
     
    // Partons in the two new systems.
    if (event[systemPtr->iParton[iFirst]].col() == event[i2New].acol()) 
      swap( i1New, i2New);
    for (int i = 0; i <= iFirst; ++i) 
      iNewSys1.push_back( systemPtr->iParton[i]);
    iNewSys1.push_back( i1New);
    iNewSys2.push_back( i2New);
    for (int i = iSecond; i < int(systemPtr->size()); ++i) 
      iNewSys2.push_back( systemPtr->iParton[i]);

  // If several gluons between sparticles then find lowest-mass gluon pair
  // and replace it by a q-qbar pair.
  } else {

    // Find lowest-mass gluon pair and adjust effective quark mass.
    int    iMin  = 0;
    int    i1Old = 0;
    int    i2Old = 0;
    double mMin  = 1e20;
    for (int i = iFirst + 1; i < iSecond - 1; ++i) { 
      int    i1Tmp = systemPtr->iParton[i];
      int    i2Tmp = systemPtr->iParton[i + 1];
      double mTmp  = (event[i1Tmp].p() + event[i2Tmp].p()).mCalc();
      if (mTmp < mMin) {
        iMin  = i;
        i1Old = i1Tmp;
        i2Old = i2Tmp;
        mMin  = mTmp;
      }
    }
    double mEff = min( mNewQ + mOffsetCloudRH, 0.4 * mMin);

    // New kinematics  by sharing mHat appropriately.
    Vec4 p1New, p2New;
    if ( !newKin( event[i1Old].p(), event[i2Old].p(), 
      mEff, mEff, p1New, p2New, false) ) return false; 

    // Insert new partons and update history.
    int i1New, i2New;
    if (event[systemPtr->iParton[0]].acol() == 0) {
      i1New = event.append( -idNewQ, 101, i1Old, 0, 0, 0, 
        0, event[i1Old].acol(), p1New, mEff, 0.);
      i2New = event.append(  idNewQ, 101, i2Old, 0, 0, 0, 
        event[i2Old].col(), 0, p2New, mEff, 0.);
    } else {
      i1New = event.append(  idNewQ, 101, i1Old, 0, 0, 0, 
        event[i1Old].col(), 0, p1New, mEff, 0.);
      i2New = event.append( -idNewQ, 101, i2Old, 0, 0, 0, 
        0, event[i2Old].acol(), p2New, mEff, 0.);
    } 
    event[i1Old].statusNeg();
    event[i2Old].statusNeg();
    event[i1Old].daughters( i1New, 0);
    event[i2Old].daughters( i2New, 0);
     
    // Partons in the two new systems.
    for (int i = 0; i < iMin; ++i) 
      iNewSys1.push_back( systemPtr->iParton[i]);
    iNewSys1.push_back( i1New);
    iNewSys2.push_back( i2New);
    for (int i = iMin + 2; i < int(systemPtr->size()); ++i) 
      iNewSys2.push_back( systemPtr->iParton[i]);
  }

  // Erase the old system and insert the two new ones instead.
  colConfig.erase(iSys);
  colConfig.insert( iNewSys1, event);
  colConfig.insert( iNewSys2, event);   

  // Done. 
  return true;

}

//--------------------------------------------------------------------------

// Produce a R-hadron from a squark and another string end.

bool RHadrons::produceSquark( ColConfig& colConfig, Event& event) {

  // Initial values.
  int    nBody  = 0;
  int    iRNow  = 0;
  int    iNewQ  = 0;
  int    iNewL  = 0;

  // Check at which end of the string the squark is located.
  int    idAbsTop = event[ systemPtr->iParton[0] ].idAbs(); 
  bool   sqAtTop  = (allowRSb && idAbsTop == idRSb) 
                 || (allowRSt && idAbsTop == idRSt);

  // Copy down system consecutively from squark end.
  int    iBeg = event.size();
  iCreRHad[iRHad] = iBeg;
  if (sqAtTop) for (int i = 0; i < systemPtr->size(); ++i) 
    event.copy( systemPtr->iParton[i], 102);
  else         for (int i = systemPtr->size() - 1; i >= 0; --i) 
    event.copy( systemPtr->iParton[i], 102);
  int    iEnd = event.size() - 1; 

  // Input flavours. From now on H = heavy and L = light string ends. 
  int    idOldH = event[iBeg].id(); 
  int    idOldL = event[iEnd].id();

  // Pick new flavour to form R-hadron. 
  FlavContainer flavOld( idOldH%10);
  int    idNewQ = flavSelPtr->pick(flavOld).id;
  int    idRHad = toIdWithSquark( idOldH, idNewQ);
  if (idRHad == 0) {
     infoPtr->errorMsg("Error in RHadrons::produceSquark: "
       "cannot form R-hadron code");
     return false;
  }

  // Target mass of R-hadron and z value of fragmentation function.
  double mRHad  = particleDataPtr->m0(idRHad) + event[iBeg].m() 
                - ( (abs(idOldH) == idRSb) ? m0Sb : m0St );
  double z      = zSelPtr->zFrag( idOldH, idNewQ, mRHad*mRHad);

  // Basic kinematics of string piece where break is to occur.
  Vec4   pOldH  = event[iBeg].p();
  int    iOldL  = iBeg + 1;
  Vec4   pOldL  = event[iOldL].p();
  double mOldL  = event[iOldL].m();
  double mNewH  = mRHad / z;
  double sSys   = (pOldH + pOldL).m2Calc();
  double sRem   = (1. - z) * (sSys - mNewH*mNewH);
  double sMin   = pow2(mOldL + mCollapseRH); 

  // If too low remaining mass in system then add one more parton to it. 
  while ( ( sRem < sMin || sSys < pow2(mNewH + mOldL + MSAFETY) )
    && iOldL < iEnd ) {    
    ++iOldL;
    pOldL      += event[iOldL].p();
    mOldL       = event[iOldL].m();
    sSys        = (pOldH + pOldL).m2Calc();
    sRem        = (1. - z) * (sSys - mNewH*mNewH);
    sMin        = pow2(mOldL + mCollapseRH); 
  }

  // If enough mass then split off R-hadron and reduced system.
  if ( sRem > sMin && sSys > pow2(mNewH + mOldL + MSAFETY) ) {
    Vec4 pNewH, pNewL;
    if ( !newKin( pOldH, pOldL, mNewH, mOldL, pNewH, pNewL) ) {
      infoPtr->errorMsg("Error in RHadrons::produceSquark: "
       "failed to construct kinematics with reduced system");
      return false;
    }

    // Insert R-hadron with new momentum. 
    iRNow       = event.append( idRHad, 104, iBeg, iOldL, 0, 0, 0, 0,
      z * pNewH, mRHad, 0.);
 
    // Reduced system with new string endpoint and modified recoiler. 
    idNewQ      = -idNewQ;
    bool hasCol = (idNewQ > 0 && idNewQ < 10) || idNewQ < -10;
    int  col    = (hasCol) ? event[iOldL].acol() : 0;
    int  acol   = (hasCol) ? 0 : event[iOldL].col(); 
    iNewQ       = event.append( idNewQ, 105, iBeg, iOldL, 0, 0, col, acol,
      (1. - z) * pNewH, (1. - z) * mNewH, 0.);
    iNewL       = event.copy( iOldL, 105);
    event[iNewL].mothers( iBeg, iOldL);
    event[iNewL].p( pNewL);

    // Done with processing of split to R-hadron and reduced system.
    nBody = 3;
  }

  // Two-body final state: form light hadron from endpoint and new flavour.
  if (nBody == 0) {
    FlavContainer flav1( idOldL);
    FlavContainer flav2( -idNewQ);
    int iTry   = 0;
    int idNewL = flavSelPtr->combine( flav1, flav2);
    while (++iTry < NTRYMAX && idNewL == 0) 
      idNewL = flavSelPtr->combine( flav1, flav2);    
    if (idNewL == 0) {
       infoPtr->errorMsg("Error in RHadrons::produceSquark: "
         "cannot form light hadron code");
       return false;
    }
    double mNewL = particleDataPtr->mass( idNewL); 
    
    // Check that energy enough for light hadron and R-hadron.
    if ( sSys > pow2(mRHad + mNewL + MSAFETY) ) { 
      Vec4 pRHad, pNewL;
      if ( !newKin( pOldH, pOldL, mRHad, mNewL, pRHad, pNewL) ) {
        infoPtr->errorMsg("Error in RHadrons::produceSquark: "
         "failed to construct kinematics for two-hadron decay");
        return false;
      }

      // Insert R-hadron and light hadron. 
      iRNow = event.append( idRHad, 104, iBeg, iOldL, 0, 0, 0, 0,
        pRHad, mRHad, 0.);
      event.append( idNewL, 105, iBeg, iOldL, 0, 0, 0, 0, 
        pNewL, mNewL, 0.);
   
      // Done for two-body case.
      nBody = 2;
    }
  }

  // Final case: for very low mass collapse to one single R-hadron.  
  if (nBody == 0) { 
    idRHad = toIdWithSquark( idOldH, idOldL);
    if (idRHad == 0) {
       infoPtr->errorMsg("Error in RHadrons::produceSquark: "
         "cannot form R-hadron code");
       return false;
    }

    // Insert R-hadron with new momentum. 
    iRNow = event.append( idRHad, 104, iBeg, iOldL, 0, 0, 0, 0,
      systemPtr->pSum, systemPtr->mass, 0.);

    // Done with one-body case.
    nBody = 1;
  }
      
  // History bookkeeping: the R-hadron and processed partons. 
  iRHadron[iRHad] = iRNow;
  int iLast = event.size() - 1;
  for (int i = iBeg; i <= iOldL; ++i) {
    event[i].statusNeg(); 
    event[i].daughters( iRNow, iLast);
  }  

  // Remove old string system and, if needed, insert a new one.
  colConfig.erase(iSys);
  if (nBody == 3) {
    vector<int> iNewSys;
    iNewSys.push_back( iNewQ);
    iNewSys.push_back( iNewL);
    for ( int i = iOldL + 1; i <= iEnd; ++i) iNewSys.push_back( i);
    colConfig.insert( iNewSys, event);
  }     

  // Copy lifetime and vertex from sparticle to R-hadron.
  event[iRNow].tau( event[iBef].tau() );
  if (event[iBef].hasVertex()) event[iRNow].vProd( event[iBef].vProd() );
 
  // Done with production of a R-hadron from a squark.  
  return true;

}

//--------------------------------------------------------------------------

// Produce a R-hadron from a gluino and two string ends (or a gluon).

bool RHadrons::produceGluino( ColConfig& colConfig, Event& event) {

  // Initial values.
  int    iGlui   = 0; 
  int    idSave  = 0; 
  int    idQLeap = 0;
  bool   isDiq1  = false;
  vector<int> iSide1, iSide2, iSideTmp, iNewSys1, iNewSys2;
  Vec4   pGlui, pSide1, pSide2;

  // Decide whether to produce a gluinoball or not.
  bool isGBall = (rndmPtr->flat() < probGluinoballRH);
         
  // Extract one string system on either side of the gluino.
  for (int i = 0; i < systemPtr->size(); ++i) {
    int  iTmp   = systemPtr->iParton[i];
    if (iGlui == 0 && event[ iTmp ].id() == idRGo) {
      iGlui     = iTmp;
      pGlui     = event[ iTmp ].p();
    } else if (iGlui == 0) {
      iSideTmp.push_back( iTmp);
      pSide1   += event[ iTmp ].p();
    } else {
      iSide2.push_back( iTmp);
      pSide2   += event[ iTmp ].p();
    }
  }
      
  // Order sides from gluino outwards and with lowest relative mass first.
  for (int i = int(iSideTmp.size()) - 1; i >= 0; --i) 
    iSide1.push_back( iSideTmp[i]);
  double m1H    = (pGlui + pSide1).mCalc() - event[ iSide1.back() ].m();
  double m2H    = (pGlui + pSide2).mCalc() - event[ iSide2.back() ].m();
  if (m2H < m1H) {
    swap( iSide1, iSide2);
    swap( pSide1, pSide2);
  }

  // Begin loop over two sides of gluino, with lowest relative mass first.
  for (int iSide = 1; iSide < 3; ++iSide) {

    // Begin processing of lower-mass string side.
    int    idRHad = 0;
    double mRHad  = 0.;
    int    nBody  = 0;
    int    iRNow  = 0;
    int    iNewQ  = 0;
    int    iNewL  = 0;
    int    statusRHad = 0;

    // Copy down system consecutively from gluino end.
    int    iBeg   = event.size();
    if (iSide == 1) {
      iCreRHad[iRHad] = iBeg;
      event.copy( iGlui, 102);
      for (int i = 0; i < int(iSide1.size()); ++i) 
        event.copy( iSide1[i], 102);
    } else {
      event.copy( iGlui, 102);
      for (int i = 0; i < int(iSide2.size()); ++i) 
        event.copy( iSide2[i], 102);
    }
    int    iEnd   = event.size() - 1; 

    // Pick new flavour to help form R-hadron. Initial colour values.
    int    idOldL = event[iEnd].id();
    int    idNewQ = 21;
    if (!isGBall) {
      do {
        FlavContainer flavOld( idOldL);
        idNewQ = -flavSelPtr->pick(flavOld).id;
      } while (iSide == 2 && isDiq1 && abs(idNewQ) > 10);
      if (iSide == 1) isDiq1 = (abs(idNewQ) > 10);
    }
    bool   hasCol = (event[iEnd].col() > 0);
    int    colR   = 0;
    int    acolR  = 0;

    // Target intermediary mass or R-hadron mass.
    if (iSide == 1) {
      idSave      = idNewQ;
      idRHad      = (hasCol) ? 1009002 : -1009002;
      if (hasCol) colR  = event[iBeg].col();
      else        acolR = event[iBeg].acol();
      statusRHad  = 103;
      double mNewQ = particleDataPtr->constituentMass( idNewQ);
      if (isGBall) mNewQ *= 0.5;
      mRHad       = event[iBeg].m() + mOffsetCloudRH + mNewQ;
    } else {
      idRHad      = toIdWithGluino( idSave, idNewQ);
      if (idRHad == 0) {
         infoPtr->errorMsg("Error in RHadrons::produceGluino: "
           "cannot form R-hadron code");
         return false;
      }
      statusRHad  = 104;
      mRHad       = particleDataPtr->m0(idRHad) + event[iBeg].m() - m0Go;
    }
      
    // z value of fragmentation function.
    double z      = zSelPtr->zFrag( idRGo, idNewQ, mRHad*mRHad);

    // Basic kinematics of string piece where break is to occur.
    Vec4   pOldH  = event[iBeg].p();
    int    iOldL  = iBeg + 1;
    Vec4   pOldL  = event[iOldL].p();
    double mOldL  = event[iOldL].m();
    double mNewH  = mRHad / z;
    double sSys   = (pOldH + pOldL).m2Calc();
    double sRem   = (1. - z) * (sSys - mNewH*mNewH);
    double sMin   = pow2(mOldL + mCollapseRH); 

    // If too low remaining mass in system then add one more parton to it. 
    while ( ( sRem < sMin || sSys < pow2(mNewH + mOldL + MSAFETY) )
      && iOldL < iEnd ) {    
      ++iOldL;
      pOldL      += event[iOldL].p();
      mOldL       = event[iOldL].m();
      sSys        = (pOldH + pOldL).m2Calc();
      sRem        = (1. - z) * (sSys - mNewH*mNewH);
      sMin        = pow2(mOldL + mCollapseRH); 
    }

    // If enough mass then split off R-hadron and reduced system.
    if ( sRem > sMin && sSys > pow2(mNewH + mOldL + MSAFETY) ) {
      Vec4 pNewH, pNewL;
      if ( !newKin( pOldH, pOldL, mNewH, mOldL, pNewH, pNewL) ) {
        infoPtr->errorMsg("Error in RHadrons::produceGluino: "
         "failed to construct kinematics with reduced system");
        return false;
      }

      // Insert intermediary or R-hadron with new momentum, less colour.
      iRNow       = event.append( idRHad, statusRHad, iBeg, iOldL, 
        0, 0, colR, acolR, z * pNewH, mRHad, 0.);
 
      // Reduced system with new string endpoint and modified recoiler. 
      if (!isGBall) idNewQ = -idNewQ;
      int  colN   = (hasCol) ? 0 : event[iOldL].acol();
      int  acolN  = (hasCol) ? event[iOldL].col() : 0; 
      iNewQ       = event.append( idNewQ, 105, iBeg, iOldL, 0, 0, 
        colN, acolN, (1. - z) * pNewH, (1. - z) * mNewH, 0.);
      iNewL       = event.copy( iOldL, 105);
      event[iNewL].mothers( iBeg, iOldL);
      event[iNewL].p( pNewL);

      // Collect partons of new string system.
      if (iSide == 1) {
        iNewSys1.push_back( iNewQ);
        iNewSys1.push_back( iNewL);
        for ( int i = iOldL + 1; i <= iEnd; ++i) iNewSys1.push_back( i);
      } else {
        iNewSys2.push_back( iNewQ);
        iNewSys2.push_back( iNewL);
        for ( int i = iOldL + 1; i <= iEnd; ++i) iNewSys2.push_back( i);
      }     

      // Done with processing of split to R-hadron and reduced system.
      nBody = 3;
    }

    // For side-1 low-mass glueball system reabsorb full momentum. 
    if (nBody == 0 && isGBall && iSide == 1) { 
      idQLeap    = event[iOldL].id();
      Vec4 pNewH = event[iBeg].p() + pOldL;

      // Insert intermediary R-hadron with new momentum, less colour.
      iRNow      = event.append( idRHad, statusRHad, iBeg, iEnd, 
        0, 0, colR, acolR, pNewH, pNewH.mCalc(), 0.);
      nBody = 1;
    }
      
    // Two-body final state: form light hadron from endpoint and new flavour.
    // Also for side 2 if gluinoball formation gives quark from side 1.
    if (nBody == 0 && (!isGBall || (iSide == 2 && idQLeap != 0) )) {
      if (isGBall) idNewQ = -idQLeap;
      FlavContainer flav1( idOldL);
      FlavContainer flav2( -idNewQ);
      int iTry   = 0;
      int idNewL = flavSelPtr->combine( flav1, flav2);
      while (++iTry < NTRYMAX && idNewL == 0) 
        idNewL = flavSelPtr->combine( flav1, flav2);    
      if (idNewL == 0) {
         infoPtr->errorMsg("Error in RHadrons::produceGluino: "
           "cannot form light hadron code");
         return false;
      }
      double mNewL = particleDataPtr->mass( idNewL); 
    
      // Check that energy enough for light hadron and R-hadron.
      if ( sSys > pow2(mRHad + mNewL + MSAFETY) ) { 
        Vec4 pRHad, pNewL;
        if ( !newKin( pOldH, pOldL, mRHad, mNewL, pRHad, pNewL) ) {
          infoPtr->errorMsg("Error in RHadrons::produceGluino: "
           "failed to construct kinematics for two-hadron decay");
          return false;
        }

        // Insert intermediary or R-hadron and light hadron. 
        iRNow = event.append( idRHad, statusRHad, iBeg, iOldL, 0, 0,
          colR, acolR, pRHad, mRHad, 0.);
        event.append( idNewL, 105, iBeg, iOldL, 0, 0, 0, 0, 
          pNewL, mNewL, 0.);
   
        // Done for two-body case.
        nBody   = 2;
        // For this case gluinoball should be handled as normal flavour.
        isGBall = false;
      }
    }

    // Final case: for very low mass collapse to one single R-hadron.  
    if (nBody == 0 && (!isGBall || (iSide == 2 && idQLeap != 0) )) { 
      if (isGBall) idSave = idQLeap;
      if (iSide == 1) idSave = idOldL;
      else            idRHad = toIdWithGluino( idSave, idOldL);
      if (idRHad == 0) {
         infoPtr->errorMsg("Error in RHadrons::produceGluino: "
           "cannot form R-hadron code");
         return false;
      }

      // Insert R-hadron with new momentum. 
      iRNow = event.append( idRHad, statusRHad, iBeg, iOldL, 0, 0, 
        colR, acolR, pOldH + pOldL, (pOldH + pOldL).mCalc(), 0.);

      // Done with one-body case.
      // Even if hoped-for, it was not possible to create a gluinoball.
      isGBall = false;
    }
      
    // History bookkeeping: the processed partons. 
    int iLast = event.size() - 1;
    for (int i = iBeg; i <= iOldL; ++i) {
      event[i].statusNeg(); 
      event[i].daughters( iRNow, iLast);
    }  

    // End of loop over two sides of the gluino.
    iGlui   = iRNow;
  }

  // History bookkeeping: insert R-hadron; delete old string system. 
  if (iGlui == 0) {
     infoPtr->errorMsg("Error in RHadrons::produceGluino: "
           "could not handle gluinoball kinematics");
     return false;
  }
  iRHadron[iRHad] = iGlui;
  colConfig.erase(iSys);

  // Non-gluinoball: insert (up to) two new string systems.
  if (!isGBall) {
    if (iNewSys1.size() > 0) colConfig.insert( iNewSys1, event);
    if (iNewSys2.size() > 0) colConfig.insert( iNewSys2, event);

  // Gluinoball with enough energy in first string: 
  // join two temporary gluons into one. 
  } else if (idQLeap == 0) {
    int iG1   = iNewSys1[0];
    int iG2   = iNewSys2[0];
    int colG  = event[iG1].col()  + event[iG2].col();  
    int acolG = event[iG1].acol() + event[iG2].acol();  
    Vec4 pG   = event[iG1].p()    + event[iG2].p(); 
    int iG12  = event.append( 21, 107, iG1, iG2, 0, 0, colG, acolG, 
      pG, pG.mCalc(), 0.);

    // Temporary gluons no longer needed, but new colour to avoid warnings.
    event[iG1].id( 21);
    event[iG2].id( 21);
    event[iG1].statusNeg();
    event[iG2].statusNeg();
    event[iG1].daughter1( iG12);
    event[iG2].daughter1( iG12);
    int colBridge = event.nextColTag();
    if (event[iG1].col() == 0) {
      event[iG1].col(  colBridge);
      event[iG2].acol( colBridge);
    } else {
      event[iG1].acol( colBridge);
      event[iG2].col(  colBridge);
    } 

    // Form the remnant system from which the R-hadron has been split off. 
    vector<int> iNewSys12;
    for (int i = int(iNewSys1.size()) - 1; i > 0; --i) 
      iNewSys12.push_back( iNewSys1[i]);
    iNewSys12.push_back( iG12);
    for (int i = 1; i < int(iNewSys2.size()); ++i) 
      iNewSys12.push_back( iNewSys2[i]);
    colConfig.insert( iNewSys12, event); 

  // Gluinoball where side 1 was fully eaten, and its flavour content
  // should leap over to the other side, to a gluon there.
  } else {
    int iG2   = iNewSys2[0];
    event[iG2].id( idQLeap);
    colConfig.insert( iNewSys2, event);
  }

  // Copy lifetime and vertex from sparticle to R-hadron.
  event[iGlui].tau( event[iBef].tau() );
  if (event[iBef].hasVertex()) event[iGlui].vProd( event[iBef].vProd() );
 
  // Done with production of a R-hadron from a gluino.  
  return true;

}

//--------------------------------------------------------------------------

// Form a R-hadron code from a squark and a (di)quark code.
// First argument is the (anti)squark, second the (anti)(di)quark.

int RHadrons::toIdWithSquark( int id1, int id2) {

  // Check that physical combination; return 0 if not.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  if (id2Abs < 10 && id1 > 0 && id2 > 0) return 0;
  if (id2Abs < 10 && id1 < 0 && id2 < 0) return 0;
  if (id2Abs > 10 && id1 > 0 && id2 < 0) return 0;
  if (id2Abs > 10 && id1 < 0 && id2 > 0) return 0;

  // Form R-hadron code. Flip sign for antisquark. 
  bool isSt = (id1Abs == idRSt);
  int idRHad = 1000000;
  if (id2Abs < 10) idRHad += ((isSt) ? 600 : 500) + 10 * id2Abs + 2;
  else idRHad += ((isSt) ? 6000 : 5000) + 10 * (id2Abs/100) + id2Abs%10;
  if (id1 < 0) idRHad = -idRHad;

  // Done.
  return idRHad;

}

//--------------------------------------------------------------------------

// Resolve a R-hadron code into a squark and a (di)quark.
// Return a pair, where first is the squark and the second the (di)quark.

pair<int,int> RHadrons::fromIdWithSquark( int idRHad) {

  // Find squark flavour content. 
  int idLight = (abs(idRHad) - 1000000) / 10;
  int idSq    = (idLight < 100) ? idLight/10 : idLight/100;
  int id1     = (idSq == 6) ? idRSt : idRSb;
  if (idRHad < 0) id1 = -id1;

  // Find light (di)quark flavour content. 
  int id2     =  (idLight < 100) ? idLight%10 : idLight%100;
  if (id2 > 10) id2 = 100 * id2 + abs(idRHad)%10;
  if ((id2 < 10 && idRHad > 0) || (id2 > 10 && idRHad < 0)) id2 = -id2;

  // Done.
  return make_pair( id1, id2);

}   
 
//--------------------------------------------------------------------------

// Form a R-hadron code from two quark/diquark endpoints and a gluino.

int RHadrons::toIdWithGluino( int id1, int id2) {

  // Check that physical combination; return 0 if not. Handle gluinoball.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  if (id1Abs == 21 && id2Abs == 21) return 1000993;
  int idMax  = max( id1Abs, id2Abs);
  int idMin  = min( id1Abs, id2Abs);
  if (idMin > 10) return 0;
  if (idMax > 10 && id1 > 0 && id2 < 0) return 0;
  if (idMax > 10 && id1 < 0 && id2 > 0) return 0;
  if (idMax < 10 && id1 > 0 && id2 > 0) return 0;
  if (idMax < 10 && id1 < 0 && id2 < 0) return 0;

  // Form R-meson code. Flip sign for antiparticle.
  int idRHad = 0;
  if (idMax < 10) {
    idRHad = 1009003 + 100 * idMax + 10 * idMin;
    if (idMin != idMax && idMax%2 == 1) {
      if (id1Abs == idMax && id1 > 0) idRHad = -idRHad;
      if (id2Abs == idMax && id2 > 0) idRHad = -idRHad;
    }
    if (idMin != idMax && idMax%2 == 0) {
      if (id1Abs == idMax && id1 < 0) idRHad = -idRHad;
      if (id2Abs == idMax && id2 < 0) idRHad = -idRHad;
    }

  // Form R-baryon code. Flip sign for antiparticle.
  } else {
    int idA = idMax/1000;
    int idB = (idMax/100)%10;
    int idC = idMin;
    if (idC > idB) swap( idB, idC);
    if (idB > idA) swap( idA, idB);    
    if (idC > idB) swap( idB, idC);
    idRHad  = 1090004 + 1000 * idA + 100 * idB + 10 * idC;
    if (id1 < 0) idRHad = -idRHad;
  }

  // Done.
  return idRHad;

}

//--------------------------------------------------------------------------

// Resolve a R-hadron code into two quark/diquark endpoints (and a gluino).
// Return a pair, where first carries colour and second anticolour.

pair<int,int> RHadrons::fromIdWithGluino( int idRHad) {

  // Find light flavour content of R-hadron.
  int idLight = (abs(idRHad) - 1000000) / 10; 
  int id1, id2, idTmp, idA, idB, idC; 

  // Gluinoballs: split g into d dbar or u ubar.
  if (idLight < 100) {
    id1 = (rndmPtr->flat() < 0.5) ? 1 : 2;
    id2 = -id1;

  // Gluino-meson: split into q + qbar.
  } else if (idLight < 1000) {
    id1 = (idLight / 10) % 10;  
    id2 = -(idLight % 10);
    // Flip signs when first quark of down-type.
    if (id1%2 == 1) {
      idTmp = id1;
      id1   = -id2;
      id2   = -idTmp;
    }

  // Gluino-baryon: split to q + qq (diquark). 
  // Pick diquark at random, except if c or b involved.
  } else {
    idA = (idLight / 100) % 10;
    idB = (idLight / 10) % 10;
    idC = idLight % 10;
    double rndmQ = 3. * rndmPtr->flat();
    if (idA > 3) rndmQ = 0.5;
    if (rndmQ < 1.) {
      id1 = idA;
      id2 = 1000 * idB + 100 * idC + 3;
      if (idB != idC && rndmPtr->flat() > diquarkSpin1RH) id2 -= 2; 
    } else if (rndmQ < 2.) {
      id1 = idB;
      id2 = 1000 * idA + 100 * idC + 3;
      if (idA != idC && rndmPtr->flat() > diquarkSpin1RH) id2 -= 2; 
    } else {
      id1 = idC;
      id2 = 1000 * idA + 100 * idB +3;
      if (idA != idB && rndmPtr->flat() > diquarkSpin1RH) id2 -= 2; 
    }
  }

  // Flip signs for anti-R-hadron.
  if (idRHad < 0) {
    idTmp = id1;
    id1   = -id2;
    id2   = -idTmp;
  }

  // Done.
  return make_pair( id1, id2);

}   
 
//--------------------------------------------------------------------------

// Construct modified four-vectors to match modified masses:
// minimal reshuffling of momentum along common axis. 
// Note that last two arguments are used to provide output!

bool RHadrons::newKin( Vec4 pOld1, Vec4 pOld2, double mNew1, double mNew2,
  Vec4& pNew1, Vec4& pNew2, bool checkMargin) {

  // Squared masses in initial and final kinematics.
  double sSum   = (pOld1 + pOld2).m2Calc();
  double sOld1  = pOld1.m2Calc();
  double sOld2  = pOld2.m2Calc();
  double sNew1  = mNew1 * mNew1;
  double sNew2  = mNew2 * mNew2;

  // Check that kinematically possible.
  if (checkMargin && pow2(mNew1 + mNew2 + MSAFETY) > sSum) return false;

  // Transfer coefficients to give four-vectors with the new masses.
  double lamOld = sqrt( pow2(sSum - sOld1 - sOld2) - 4. * sOld1 * sOld2 );
  double lamNew = sqrt( pow2(sSum - sNew1 - sNew2) - 4. * sNew1 * sNew2 );   
  double move1  = (lamNew * (sSum - sOld1 + sOld2) 
                -  lamOld * (sSum - sNew1 + sNew2)) / (2. * sSum * lamOld);
  double move2  = (lamNew * (sSum + sOld1 - sOld2) 
                -  lamOld * (sSum + sNew1 - sNew2)) / (2. * sSum * lamOld);
  
  // Construct final vectors. Done.
  pNew1 = (1. + move1) * pOld1 - move2 * pOld2;
  pNew2 = (1. + move2) * pOld2 - move1 * pOld1;
  return true;

}   
 
//==========================================================================

} // end namespace Pythia8
