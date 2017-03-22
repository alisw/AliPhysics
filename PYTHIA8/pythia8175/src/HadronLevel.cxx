// HadronLevel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the HadronLevel class.

#include "HadronLevel.h"

namespace Pythia8 {

//==========================================================================

// The HadronLevel class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// For breaking J-J string, pick a Gamma by taking a step with fictitious mass.
const double HadronLevel::JJSTRINGM2MAX  = 25.; 
const double HadronLevel::JJSTRINGM2FRAC = 0.1; 

// Iterate junction rest frame boost until convergence or too many tries.
const double HadronLevel::CONVJNREST     = 1e-5;
const int HadronLevel::NTRYJNREST        = 20; 

// Typical average transvere primary hadron mass <mThad>. 
const double HadronLevel::MTHAD          = 0.9; 

//--------------------------------------------------------------------------

// Find settings. Initialize HadronLevel classes as required.

bool HadronLevel::init(Info* infoPtrIn, Settings& settings, 
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn, 
  Couplings* couplingsPtrIn, TimeShower* timesDecPtr, 
  RHadrons* rHadronsPtrIn, DecayHandler* decayHandlePtr, 
  vector<int> handledParticles) {

  // Save pointers.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  couplingsPtr    = couplingsPtrIn;
  rHadronsPtr     = rHadronsPtrIn; 

  // Main flags.
  doHadronize     = settings.flag("HadronLevel:Hadronize");
  doDecay         = settings.flag("HadronLevel:Decay");
  doBoseEinstein  = settings.flag("HadronLevel:BoseEinstein");

  // Boundary mass between string and ministring handling.
  mStringMin      = settings.parm("HadronLevel:mStringMin");

  // For junction processing.
  eNormJunction   = settings.parm("StringFragmentation:eNormJunction");

  // Allow R-hadron formation.
  allowRH         = settings.flag("RHadrons:allow");

  // Particles that should decay or not before Bose-Einstein stage.
  widthSepBE      = settings.parm("BoseEinstein:widthSep");

  // Hadron scattering --rjc
  doHadronScatter = settings.flag("HadronScatter:scatter");
  hsAfterDecay    = settings.flag("HadronScatter:afterDecay");

  // Initialize auxiliary fragmentation classes.
  flavSel.init(settings, rndmPtr);
  pTSel.init(settings, *particleDataPtr, rndmPtr);
  zSel.init(settings, *particleDataPtr, rndmPtr);

  // Initialize auxiliary administrative class.
  colConfig.init(infoPtr, settings, &flavSel);

  // Initialize string and ministring fragmentation.
  stringFrag.init(infoPtr, settings, particleDataPtr, rndmPtr, 
    &flavSel, &pTSel, &zSel);
  ministringFrag.init(infoPtr, settings, particleDataPtr, rndmPtr, 
    &flavSel, &pTSel, &zSel);
 
  // Initialize particle decays.  
  decays.init(infoPtr, settings, particleDataPtr, rndmPtr, couplingsPtr, 
    timesDecPtr, &flavSel, decayHandlePtr, handledParticles); 

  // Initialize BoseEinstein. 
  boseEinstein.init(infoPtr, settings, *particleDataPtr); 

  // Initialize HadronScatter --rjc
  if (doHadronScatter)
    hadronScatter.init(infoPtr, settings, rndmPtr, particleDataPtr);

  // Initialize Hidden-Valley fragmentation, if necessary.
  useHiddenValley = hiddenvalleyFrag.init(infoPtr, settings, 
    particleDataPtr, rndmPtr);

  // Send flavour and z selection pointers to R-hadron machinery.
  rHadronsPtr->fragPtrs( &flavSel, &zSel);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Hadronize and decay the next parton-level.

bool HadronLevel::next( Event& event) {

  // Do Hidden-Valley fragmentation, if necessary.
  if (useHiddenValley) hiddenvalleyFrag.fragment(event);

  // Colour-octet onia states must be decayed to singlet + gluon.
  if (!decayOctetOnia(event)) return false;

  // Possibility of hadronization inside decay, but then no BE second time.
  // Hadron scattering, first pass only --rjc
  bool moreToDo, firstPass = true;
  bool doBoseEinsteinNow = doBoseEinstein;
  do {
    moreToDo = false;

    // First part: string fragmentation.   
    if (doHadronize) {

      // Find the complete colour singlet configuration of the event.
      if (!findSinglets( event)) return false;

      // Fragment off R-hadrons, if necessary. 
      if (allowRH && !rHadronsPtr->produce( colConfig, event)) 
        return false;

      // Process all colour singlet (sub)system
      for (int iSub = 0; iSub < colConfig.size(); ++iSub) {

        // Collect sequentially all partons in a colour singlet subsystem.
        colConfig.collect(iSub, event);

        // String fragmentation of each colour singlet (sub)system.  
        if ( colConfig[iSub].massExcess > mStringMin ) {
          if (!stringFrag.fragment( iSub, colConfig, event)) return false; 

        // Low-mass string treated separately. Tell if diffractive system.
        } else { 
          bool isDiff = infoPtr->isDiffractiveA() 
            || infoPtr->isDiffractiveB();
          if (!ministringFrag.fragment( iSub, colConfig, event, isDiff)) 
            return false;
        } 
      }
    }

    // Hadron scattering --rjc
    if (doHadronScatter && !hsAfterDecay && firstPass)
      hadronScatter.scatter(event);

    // Second part: sequential decays of short-lived particles (incl. K0).
    if (doDecay) {
    
      // Loop through all entries to find those that should decay.
      int iDec = 0;
      do {
        Particle& decayer = event[iDec];
        if ( decayer.isFinal() && decayer.canDecay() && decayer.mayDecay() 
          && (decayer.mWidth() > widthSepBE || decayer.idAbs() == 311) ) {
          decays.decay( iDec, event); 
          if (decays.moreToDo()) moreToDo = true;
        }
      } while (++iDec < event.size());
    }

    // Hadron scattering --rjc
    if (doHadronScatter && hsAfterDecay && firstPass)
      hadronScatter.scatter(event);

    // Third part: include Bose-Einstein effects among current particles.
    if (doBoseEinsteinNow) {
      if (!boseEinstein.shiftEvent(event)) return false;
      doBoseEinsteinNow = false;
    }
    
    // Fourth part: sequential decays also of long-lived particles.
    if (doDecay) {
    
      // Loop through all entries to find those that should decay.
      int iDec = 0;
      do {
        Particle& decayer = event[iDec];
        if ( decayer.isFinal() && decayer.canDecay() && decayer.mayDecay() ) {
          decays.decay( iDec, event); 
          if (decays.moreToDo()) moreToDo = true;
        }
      } while (++iDec < event.size());
    }

  // Normally done first time around, but sometimes not (e.g. Upsilon).
  } while (moreToDo);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Allow more decays if on/off switches changed.
// Note: does not do sequential hadronization, e.g. for Upsilon.

bool HadronLevel::moreDecays( Event& event) {

  // Colour-octet onia states must be decayed to singlet + gluon.
  if (!decayOctetOnia(event)) return false;
    
  // Loop through all entries to find those that should decay.
  int iDec = 0;
  do {
    if ( event[iDec].isFinal() && event[iDec].canDecay() 
      && event[iDec].mayDecay() ) decays.decay( iDec, event); 
  } while (++iDec < event.size());

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Decay colour-octet onium states.

bool HadronLevel::decayOctetOnia(Event& event) {

  // Onium states to be decayed.
  int idOnium[6] = { 9900443, 9900441, 9910441, 
                     9900553, 9900551, 9910551 };
  
  // Loop over particles and identify onia.
  for (int iDec = 0; iDec < event.size(); ++iDec) 
  if (event[iDec].isFinal()) {
    int id = event[iDec].id();  
    bool isOnium = false;
    for (int j = 0; j < 6; ++j) if (id == idOnium[j]) isOnium = true;
    
    // Decay any onia encountered.
    if (isOnium) { 
      if (!decays.decay( iDec, event)) return false;

      // Set colour flow by hand: gluon inherits octet-onium state.
      int iGlu = event.size() - 1;
      event[iGlu].cols( event[iDec].col(), event[iDec].acol() );
    }
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Trace colour flow in the event to form colour singlet subsystems.

bool HadronLevel::findSinglets(Event& event) {
  
  // Find a list of final partons and of all colour ends and gluons.
  iColEnd.resize(0);
  iAcolEnd.resize(0);
  iColAndAcol.resize(0);
  for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
    if (event[i].col() > 0 && event[i].acol() > 0) iColAndAcol.push_back(i);
    else if (event[i].col() > 0) iColEnd.push_back(i);
    else if (event[i].acol() > 0) iAcolEnd.push_back(i); 
  }  

  // Begin arrange the partons into separate colour singlets.
  colConfig.clear();
  iPartonJun.resize(0);
  iPartonAntiJun.resize(0);

  // Junctions: loop over them, and identify kind.
  for (int iJun = 0; iJun < event.sizeJunction(); ++iJun)     
  if (event.remainsJunction(iJun)) {
    event.remainsJunction(iJun, false);
    int kindJun = event.kindJunction(iJun);
    iParton.resize(0);

    // Loop over junction legs.
    for (int iCol = 0; iCol < 3; ++iCol) {
      int indxCol = event.colJunction(iJun, iCol);    
      iParton.push_back( -(10 + 10 * iJun + iCol) );
      // Junctions: find color ends.
      if (kindJun % 2 == 1 && !traceFromAcol(indxCol, event, iJun, iCol)) 
	return false;       
      // Antijunctions: find anticolor ends.
      if (kindJun % 2 == 0 && !traceFromCol(indxCol, event, iJun, iCol)) 
	return false;      
    }

    // Reject triple- and higher-junction systems (physics not implemented).
    int otherJun = 0;
    for (int i = 0; i < int(iParton.size()); ++i) 
    if (iParton[i] < 0 && abs(iParton[i]) / 10 != iJun + 1) {
      if (otherJun == 0) otherJun = abs(iParton[i]) / 10; 
      else if (abs(iParton[i]) / 10 != otherJun) {
        infoPtr->errorMsg("Error in HadronLevel::findSinglets: "
          "too many junction-junction connections"); 
        return false;
      }
    }

    // Keep in memory a junction hooked up with an antijunction,
    // else store found single-junction system.
    int nNeg = 0;
    for (int i = 0; i < int(iParton.size()); ++i) if (iParton[i] < 0) 
      ++nNeg; 
    if (nNeg > 3 && kindJun % 2 == 1) { 
      for (int i = 0; i < int(iParton.size()); ++i) 
        iPartonJun.push_back(iParton[i]);
    } else if (nNeg > 3 && kindJun % 2 == 0) { 
      for (int i = 0; i < int(iParton.size()); ++i) 
        iPartonAntiJun.push_back(iParton[i]);
    } else {
      // A junction may be eliminated by insert if two quarks are nearby.
      int nJunOld = event.sizeJunction(); 
      if (!colConfig.insert(iParton, event)) return false;
      if (event.sizeJunction() < nJunOld) --iJun; 
    }
  }

  // Split junction-antijunction system into two, and store those.
  // (Only one system in extreme cases, and then second empty.)
  if (iPartonJun.size() > 0 && iPartonAntiJun.size() > 0) {
    if (!splitJunctionPair(event)) return false;
    if (!colConfig.insert(iPartonJun, event)) return false;
    if (iPartonAntiJun.size() > 0) 
      if (!colConfig.insert(iPartonAntiJun, event)) return false;
  // Error if only one of junction and antijuction left here.
  } else if (iPartonJun.size() > 0 || iPartonAntiJun.size() > 0) {
    infoPtr->errorMsg("Error in HadronLevel::findSinglets: "
      "unmatched (anti)junction"); 
    return false;
  } 

  // Open strings: pick up each colour end and trace to its anticolor end.
  for (int iEnd = 0; iEnd < int(iColEnd.size()); ++iEnd) {
    iParton.resize(0);
    iParton.push_back( iColEnd[iEnd] );
    int indxCol = event[ iColEnd[iEnd] ].col();    
    if (!traceFromCol(indxCol, event)) return false;

    // Store found open string system. Analyze its properties.
    if (!colConfig.insert(iParton, event)) return false;
  }

  // Closed strings : begin at any gluon and trace until back at it. 
  while (iColAndAcol.size() > 0) {
    iParton.resize(0);
    iParton.push_back( iColAndAcol[0] );
    int indxCol = event[ iColAndAcol[0] ].col();    
    int indxAcol = event[ iColAndAcol[0] ].acol();    
    iColAndAcol[0] = iColAndAcol.back();
    iColAndAcol.pop_back();
    if (!traceInLoop(indxCol, indxAcol, event)) return false;

    // Store found closed string system. Analyze its properties.
    if (!colConfig.insert(iParton, event)) return false;
  }
    
  // Done. 
  return true;

}

//--------------------------------------------------------------------------

// Trace a colour line, from a colour to an anticolour.

bool HadronLevel::traceFromCol(int indxCol, Event& event, int iJun, 
  int iCol) {

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
    for (int i = 0; i < int(iAcolEnd.size()); ++i)      
    if (event[ iAcolEnd[i] ].acol() == indxCol) {
      iParton.push_back( iAcolEnd[i] );
      indxCol = 0;
      iAcolEnd[i] = iAcolEnd.back();
      iAcolEnd.pop_back();
      hasFound = true;
      break;
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

    // In a pinch, check list of opposite-sign junction end colours.
    // Store in iParton list as -(10 + 10 * iAntiJun + iAntiLeg).
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
    infoPtr->errorMsg("Error in HadronLevel::traceFromCol: "
      "colour tracing failed"); 
    return false;
  } 

  // Done. 
  return true;

}

//--------------------------------------------------------------------------

// Trace a colour line, from an anticolour to a colour.

bool HadronLevel::traceFromAcol(int indxCol, Event& event, int iJun,
  int iCol) {

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
    for (int i = 0; i < int(iColEnd.size()); ++i) 
    if (event[ iColEnd[i] ].col() == indxCol) {
      iParton.push_back( iColEnd[i] );
      indxCol = 0;
      iColEnd[i] = iColEnd.back();
      iColEnd.pop_back();
      hasFound = true;
      break;
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

    // In a pinch, check list of opposite-sign junction end colours.
    // Store in iParton list as -(10 + 10 * iAntiJun + iLeg).
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
    infoPtr->errorMsg("Error in HadronLevel::traceFromAcol: "
      "colour tracing failed"); 
    return false;
  }

  // Done. 
  return true;

}

//--------------------------------------------------------------------------

// Trace a colour loop, from a colour back to the anticolour of the same.

bool HadronLevel::traceInLoop(int indxCol, int indxAcol, Event& event) {
    
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
    infoPtr->errorMsg("Error in HadronLevel::traceInLoop: "
      "colour tracing failed"); 
    return false; 
  } 

  // Done. 
  return true;

}

//--------------------------------------------------------------------------

// Split junction-antijunction system into two, or simplify other way.

bool HadronLevel::splitJunctionPair(Event& event) {

  // Construct separate index arrays for the three junction legs.
  int identJun = (-iPartonJun[0])/10;
  iJunLegA.resize(0);
  iJunLegB.resize(0);
  iJunLegC.resize(0);
  int leg = -1;
  for (int i = 0; i < int(iPartonJun.size()); ++ i) {
    if ( (-iPartonJun[i])/10 == identJun) ++leg;
    if (leg == 0) iJunLegA.push_back( iPartonJun[i] );
    else if (leg == 1) iJunLegB.push_back( iPartonJun[i] );
    else iJunLegC.push_back( iPartonJun[i] );
  }

  // Construct separate index arrays for the three antijunction legs.
  int identAnti = (-iPartonAntiJun[0])/10;
  iAntiLegA.resize(0);
  iAntiLegB.resize(0);
  iAntiLegC.resize(0);
  leg = -1;
  for (int i = 0; i < int(iPartonAntiJun.size()); ++ i) {
    if ( (-iPartonAntiJun[i])/10 == identAnti) ++leg;
    if (leg == 0) iAntiLegA.push_back( iPartonAntiJun[i] );
    else if (leg == 1) iAntiLegB.push_back( iPartonAntiJun[i] );
    else iAntiLegC.push_back( iPartonAntiJun[i] );
  }

  // Find interjunction legs, i.e. between junction and antijunction.
  int nMatch = 0;
  int legJun[3], legAnti[3], nGluLeg[3];
  if (iJunLegA.back() < 0) { legJun[nMatch] = 0;
    legAnti[nMatch] = (-iJunLegA.back())%10; ++nMatch;}
  if (iJunLegB.back() < 0) { legJun[nMatch] = 1;
    legAnti[nMatch] = (-iJunLegB.back())%10; ++nMatch;}
  if (iJunLegC.back() < 0) { legJun[nMatch] = 2;
    legAnti[nMatch] = (-iJunLegC.back())%10; ++nMatch;}

  // Loop over interjunction legs.
  for (int iMatch = 0; iMatch < nMatch; ++iMatch) {
    vector<int>& iJunLeg = (legJun[iMatch] == 0) ? iJunLegA
      : ( (legJun[iMatch] == 1) ? iJunLegB : iJunLegC );
    vector<int>& iAntiLeg = (legAnti[iMatch] == 0) ? iAntiLegA
      : ( (legAnti[iMatch] == 1) ? iAntiLegB : iAntiLegC );

    // Find number of gluons on each. Do nothing for now if none.
    nGluLeg[iMatch] = iJunLeg.size() + iAntiLeg.size() - 4;
    if (nGluLeg[iMatch] == 0) continue;

    // Else pick up the gluons on the interjunction leg in order.
    iGluLeg.resize(0);
    for (int i = 1; i < int(iJunLeg.size()) - 1; ++i) 
      iGluLeg.push_back( iJunLeg[i] );
    for (int i = int(iAntiLeg.size()) - 2; i > 0; --i) 
      iGluLeg.push_back( iAntiLeg[i] );

    // Remove those gluons from the junction/antijunction leg lists.
    iJunLeg.resize(1);
    iAntiLeg.resize(1);

   // Pick a new quark at random; for simplicity no diquarks.
    int idQ = flavSel.pickLightQ();
    int colQ, acolQ;

    // If one gluon on leg, split it into a collinear q-qbar pair.
    if (iGluLeg.size() == 1) { 
    
      // Store the new q qbar pair, sharing gluon colour and momentum.
      colQ = event[ iGluLeg[0] ].col();
      acolQ = event[ iGluLeg[0] ].acol();
      Vec4 pQ = 0.5 * event[ iGluLeg[0] ].p(); 
      double mQ = 0.5 * event[ iGluLeg[0] ].m(); 
      int iQ = event.append( idQ, 75, iGluLeg[0], 0, 0, 0, colQ, 0, pQ, mQ );
      int iQbar = event.append( -idQ, 75, iGluLeg[0], 0, 0, 0, 0, acolQ, 
        pQ, mQ );

      // Mark split gluon and update junction and antijunction legs.
      event[ iGluLeg[0] ].statusNeg();
      event[ iGluLeg[0] ].daughters( iQ, iQbar);    
      iJunLeg.push_back(iQ);         
      iAntiLeg.push_back(iQbar);

    // If several gluons on the string, decide which g-g region to split up.
    } else {

      // Evaluate mass-squared for all adjacent gluon pairs.
      m2Pair.resize(0);
      double m2Sum = 0.;
      for (int i = 0; i < int(iGluLeg.size()) - 1; ++i) {
        double m2Now = 0.5 * event[ iGluLeg[i] ].p() 
          * event[ iGluLeg[i + 1] ].p();  
        m2Pair.push_back(m2Now);
        m2Sum += m2Now;
      }
   
      // Pick breakup region with probability proportional to mass-squared.
      double m2Reg = m2Sum * rndmPtr->flat();
      int iReg = -1;
      do m2Reg -= m2Pair[++iReg];
      while (m2Reg > 0. && iReg < int(iGluLeg.size()) - 1); 
      m2Reg = m2Pair[iReg];

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
      Particle& gJun = event[ iGluLeg[iReg] ]; 
      Particle& gAnti = event[ iGluLeg[iReg + 1] ]; 
      gJun.statusNeg();
      gAnti.statusNeg();
      int dau1 = event.size();
      gJun.daughters(dau1, dau1 + 3);
      gAnti.daughters(dau1, dau1 + 3);
      int mother1 = min( iGluLeg[iReg], iGluLeg[iReg + 1]); 
      int mother2 = max( iGluLeg[iReg], iGluLeg[iReg + 1]); 

      // Can keep one of old colours but need one new so unambiguous.
      colQ = gJun.acol();
      acolQ = event.nextColTag(); 

      // Store copied gluons with reduced momenta.
      int iGjun = event.append( 21, 75, mother1, mother2, 0, 0, 
        gJun.col(), gJun.acol(), (1. - 0.5 * xPos) * gJun.p(),
        (1. - 0.5 * xPos) * gJun.m());
      int iGanti = event.append( 21, 75, mother1, mother2, 0, 0, 
        acolQ, gAnti.acol(), (1. - 0.5 * xNeg) * gAnti.p(),
        (1. - 0.5 * xNeg) * gAnti.m());

      // Store the new q qbar pair with remaining momenta.
      int iQ = event.append( idQ, 75, mother1, mother2, 0, 0, 
        colQ, 0, 0.5 * xNeg * gAnti.p(), 0.5 * xNeg * gAnti.m() );
      int iQbar = event.append( -idQ, 75, mother1, mother2, 0, 0, 
        0, acolQ, 0.5 * xPos * gJun.p(), 0.5 * xPos * gJun.m() );

      // Update junction and antijunction legs with gluons and quarks. 
      for (int i = 0; i < iReg; ++i) 
        iJunLeg.push_back( iGluLeg[i] );
      iJunLeg.push_back(iGjun);
      iJunLeg.push_back(iQ);
      for (int i = int(iGluLeg.size()) - 1; i > iReg + 1; --i) 
        iAntiLeg.push_back( iGluLeg[i] );
      iAntiLeg.push_back(iGanti);
      iAntiLeg.push_back(iQbar);
    }

    // Update end colours for both g -> q qbar and g g -> g g q qbar.
    event.endColJunction(identJun - 1, legJun[iMatch], colQ);
    event.endColJunction(identAnti - 1, legAnti[iMatch], acolQ);
  } 

  // Update list of interjunction legs after splittings above.  
  int iMatchUp = 0;
  while (iMatchUp < nMatch) {
    if (nGluLeg[iMatchUp] > 0) {
      for (int i = iMatchUp; i < nMatch - 1; ++i) {
        legJun[i] = legJun[i + 1];
        legAnti[i] = legAnti[i + 1];
        nGluLeg[i] = nGluLeg[i + 1];
      } --nMatch;
    } else ++iMatchUp;
  } 

  // Should not ever have three empty interjunction legs.
  if (nMatch == 3) {
    infoPtr->errorMsg("Error in HadronLevel::splitJunctionPair: "
      "three empty junction-junction legs"); 
    return false;
  }

  // If two legs are empty, then collapse system to a single string.
  if (nMatch == 2) {
    int legJunLeft = 3 - legJun[0] - legJun[1];
    int legAntiLeft = 3 - legAnti[0] - legAnti[1];
    vector<int>& iJunLeg = (legJunLeft == 0) ? iJunLegA
      : ( (legJunLeft == 1) ? iJunLegB : iJunLegC );
    vector<int>& iAntiLeg = (legAntiLeft == 0) ? iAntiLegA
      : ( (legAntiLeft == 1) ? iAntiLegB : iAntiLegC );
    iPartonJun.resize(0);
    for (int i = int(iJunLeg.size()) - 1; i > 0; --i) 
      iPartonJun.push_back( iJunLeg[i] );
    for (int i = 1; i < int(iAntiLeg.size()); ++i) 
      iPartonJun.push_back( iAntiLeg[i] );

    // Match up the colours where the strings are joined.
    int iColJoin  = iJunLeg[1];
    int iAcolJoin = iAntiLeg[1];
    event[iAcolJoin].acol( event[iColJoin].col() ); 

    // Other string system empty. Remove junctions from their list. Done.
    iPartonAntiJun.resize(0);
    event.eraseJunction( max(identJun, identAnti) - 1);
    event.eraseJunction( min(identJun, identAnti) - 1);
    return true;
  }  

  // If one leg is empty then, depending on string length, either 
  // (a) annihilate junction and antijunction into two simple strings, or 
  // (b) split the empty leg by borrowing energy from nearby legs.
  if (nMatch == 1) {

    // Identify the two external legs of either junction.
    vector<int>& iJunLeg0 = (legJun[0] == 0) ? iJunLegB : iJunLegA;
    vector<int>& iJunLeg1 = (legJun[0] == 2) ? iJunLegB : iJunLegC;
    vector<int>& iAntiLeg0 = (legAnti[0] == 0) ? iAntiLegB : iAntiLegA;
    vector<int>& iAntiLeg1 = (legAnti[0] == 2) ? iAntiLegB : iAntiLegC;

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
    for (int iJun = 0; iJun < 2; ++iJun) {
      int offset = (iJun == 0) ? 0 : 2;
         
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
        for (int i = 0; i < 4; ++i) pInRF[i].rotbst(MtoRF);

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
      if (iJun == 0) {
        MtoJRF = MtoRF;
        for (int i = 0; i < 3; ++i) pInJRF[i] = pInRF[i];
      } else { 
        MtoARF = MtoRF;
        for (int i = 0; i < 3; ++i) pInARF[i] = pInRF[i];
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
    for (int i = 0; i < 3; ++i) { 
      pLabJ[i] = pInJRF[i];
      pLabJ[i].rotbst(MfromJRF);
      pLabA[i] = pInARF[i];
      pLabA[i].rotbst(MfromARF);
    }

    // Calculate Lambda-measure length of three possible topologies.
    double vJvA = vJun * vAnti;
    double vJvAe2y = vJvA + sqrt(vJvA*vJvA - 1.);
    double LambdaJA = (2. * pInJRF[0].e()) * (2. * pInJRF[1].e())   
      * (2. * pInARF[0].e()) * (2. * pInARF[1].e()) * vJvAe2y;
    double Lambda00 = (2. * pLabJ[0] * pLabA[0])
      * (2. * pLabJ[1] * pLabA[1]);
    double Lambda01 = (2. * pLabJ[0] * pLabA[1])
      * (2. * pLabJ[1] * pLabA[0]);

    // Case when either topology without junctions is the shorter one.
    if (LambdaJA > min( Lambda00, Lambda01)) {  
      vector<int>& iAntiMatch0 = (Lambda00 < Lambda01) 
        ? iAntiLeg0 : iAntiLeg1;
      vector<int>& iAntiMatch1 = (Lambda00 < Lambda01) 
        ? iAntiLeg1 : iAntiLeg0;

      // Define two quark-antiquark strings.
      iPartonJun.resize(0);
      for (int i = int(iJunLeg0.size()) - 1; i > 0; --i) 
        iPartonJun.push_back( iJunLeg0[i] );
      for (int i = 1; i < int(iAntiMatch0.size()); ++i) 
        iPartonJun.push_back( iAntiMatch0[i] );
      iPartonAntiJun.resize(0);
      for (int i = int(iJunLeg1.size()) - 1; i > 0; --i) 
        iPartonAntiJun.push_back( iJunLeg1[i] );
      for (int i = 1; i < int(iAntiMatch1.size()); ++i) 
        iPartonAntiJun.push_back( iAntiMatch1[i] );

      // Match up the colours where the strings are joined.
      int iColJoin  = iJunLeg0[1];
      int iAcolJoin = iAntiMatch0[1];
      event[iAcolJoin].acol( event[iColJoin].col() ); 
      iColJoin  = iJunLeg1[1];
      iAcolJoin = iAntiMatch1[1];
      event[iAcolJoin].acol( event[iColJoin].col() ); 
 
      // Remove junctions from their list. Done.
      event.eraseJunction( max(identJun, identAnti) - 1);
      event.eraseJunction( min(identJun, identAnti) - 1);
      return true;
    }

    // Case where junction and antijunction to be separated.
    // Shuffle (p+/p-)  momentum of order <mThad> between systems,
    // times 2/3 for 120 degree in JRF, times 1/2 for two legs,
    // but not more than half of what nearest parton carries.
    double eShift = MTHAD / (3. * sqrt(vJvAe2y));
    double fracJ0 = min(0.5, eShift / pInJRF[0].e());
    double fracJ1 = min(0.5, eShift / pInJRF[0].e());
    Vec4 pFromJun = fracJ0 * pJunLeg0 + fracJ1 * pJunLeg1;  
    double fracA0 = min(0.5, eShift / pInARF[0].e());
    double fracA1 = min(0.5, eShift / pInARF[0].e());
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
    event.endColJunction(identJun - 1, legJun[0], colQ);
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
    event.endColJunction(identAnti - 1, legAnti[0], acolQ);
    int iNewA = event.append( -idQ, 76, mother1, mother2, 0, 0, 
      0, acolQ, pFromJun, pFromJun.mCalc() );
    event[mother1].daughters( iNew1, iNewA);
    event[mother2].daughters( iNew1, iNewA);
    event[iNew1].mothers( mother1, mother2);    
    event[iNew2].mothers( mother1, mother2);    

    // Bookkeep new quark and antiquark on third legs.
    if (legJun[0] == 0) iJunLegA[1] = iNewJ;
    else if (legJun[0] == 1) iJunLegB[1] = iNewJ;
    else iJunLegC[1] = iNewJ;
    if (legAnti[0] == 0) iAntiLegA[1] = iNewA;
    else if (legAnti[0] == 1) iAntiLegB[1] = iNewA;
    else iAntiLegC[1] = iNewA;

  // Done with splitting junction from antijunction.
  }
   
  // Put together new junction parton list.
  iPartonJun.resize(0);
  for (int i = 0; i < int(iJunLegA.size()); ++i) 
    iPartonJun.push_back( iJunLegA[i] );
  for (int i = 0; i < int(iJunLegB.size()); ++i) 
    iPartonJun.push_back( iJunLegB[i] );
  for (int i = 0; i < int(iJunLegC.size()); ++i) 
    iPartonJun.push_back( iJunLegC[i] );

  // Put together new antijunction parton list.
  iPartonAntiJun.resize(0);
  for (int i = 0; i < int(iAntiLegA.size()); ++i) 
    iPartonAntiJun.push_back( iAntiLegA[i] );
  for (int i = 0; i < int(iAntiLegB.size()); ++i) 
    iPartonAntiJun.push_back( iAntiLegB[i] );
  for (int i = 0; i < int(iAntiLegC.size()); ++i) 
    iPartonAntiJun.push_back( iAntiLegC[i] );

  // Now the two junction systems are separated and can be stored.
  return true;

}

//==========================================================================

} // end namespace Pythia8

