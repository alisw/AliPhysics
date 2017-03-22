// SpaceShower.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// SpaceShower class.

#include "SpaceShower.h"

namespace Pythia8 {

//==========================================================================

// The SpaceShower class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Leftover companion can give PDF > 0 at small Q2 where other PDF's = 0,
// and then one can end in infinite loop of impossible kinematics.
const int    SpaceShower::MAXLOOPTINYPDF = 10; 

// Switch to alternative (but equivalent) backwards evolution for
// g -> Q Qbar (Q = c or b) when below QTHRESHOLD * mQ2.
const double SpaceShower::CTHRESHOLD     = 2.0; 
const double SpaceShower::BTHRESHOLD     = 2.0; 

// Renew evaluation of PDF's when the pT2 step is bigger than this 
// (in addition to initial scale and c and b thresholds.)
const double SpaceShower::EVALPDFSTEP    = 0.1;

// Lower limit on PDF value in order to avoid division by zero.
const double SpaceShower::TINYPDF        = 1e-10;

// Lower limit on estimated evolution rate, below which stop.
const double SpaceShower::TINYKERNELPDF  = 1e-6;

// Lower limit on pT2, below which branching is rejected. 
const double SpaceShower::TINYPT2        = 0.25e-6;

// No attempt to do backwards evolution of a heavy (c or b) quark 
// if evolution starts at a scale pT2 < HEAVYPT2EVOL * mQ2.
const double SpaceShower::HEAVYPT2EVOL   = 1.1;

// No attempt to do backwards evolution of a heavy (c or b) quark 
// if evolution starts at a  x > HEAVYXEVOL * x_max, where 
// x_max is the largest possible x value for a g -> Q Qbar branching.
const double SpaceShower::HEAVYXEVOL     = 0.9;
  
// When backwards evolution Q -> g + Q creates a heavy quark Q,
// an earlier branching g -> Q + Qbar will restrict kinematics
// to  M_{Q Qbar}^2 > EXTRASPACEQ * 4 m_Q^2. (Smarter to be found??) 
const double SpaceShower::EXTRASPACEQ    = 2.0;

// Never pick pT so low that alphaS is evaluated too close to Lambda_3. 
const double SpaceShower::LAMBDA3MARGIN  = 1.1;

// Do not warn for large PDF ratios at small pT2 scales.
const double SpaceShower::PT2MINWARN = 1.;

// Cutoff for f_e^e at x < 1 - 10^{-10} to be used in z selection.
// Note: the x_min quantity come from 1 - x_max.
const double SpaceShower::LEPTONXMIN     = 1e-10;
const double SpaceShower::LEPTONXMAX     = 1. - 1e-10;

// Stop l -> l gamma evolution slightly above m2l.
const double SpaceShower::LEPTONPT2MIN   = 1.2;

// Enhancement of l -> l gamma trial rate to compensate imperfect modelling.
const double SpaceShower::LEPTONFUDGE    = 10.;

//--------------------------------------------------------------------------

// Initialize alphaStrong, alphaEM and related pTmin parameters.

void SpaceShower::init( BeamParticle* beamAPtrIn, 
  BeamParticle* beamBPtrIn) {

  // Store input pointers for future use. 
  beamAPtr        = beamAPtrIn;
  beamBPtr        = beamBPtrIn;

  // Main flags to switch on and off branchings.
  doQCDshower     = settingsPtr->flag("SpaceShower:QCDshower");
  doQEDshowerByQ  = settingsPtr->flag("SpaceShower:QEDshowerByQ");
  doQEDshowerByL  = settingsPtr->flag("SpaceShower:QEDshowerByL");

  // Matching in pT of hard interaction to shower evolution.
  pTmaxMatch      = settingsPtr->mode("SpaceShower:pTmaxMatch"); 
  pTdampMatch     = settingsPtr->mode("SpaceShower:pTdampMatch"); 
  pTmaxFudge      = settingsPtr->parm("SpaceShower:pTmaxFudge"); 
  pTmaxFudgeMPI   = settingsPtr->parm("SpaceShower:pTmaxFudgeMPI"); 
  pTdampFudge     = settingsPtr->parm("SpaceShower:pTdampFudge"); 

  // Optionally force emissions to be ordered in rapidity/angle.
  doRapidityOrder = settingsPtr->flag("SpaceShower:rapidityOrder");

  // Charm, bottom and lepton mass thresholds.
  mc              = particleDataPtr->m0(4); 
  mb              = particleDataPtr->m0(5); 
  m2c             = pow2(mc);
  m2b             = pow2(mb);

  // Parameters of scale choices.
  renormMultFac     = settingsPtr->parm("SpaceShower:renormMultFac");
  factorMultFac     = settingsPtr->parm("SpaceShower:factorMultFac");

  // Parameters of alphaStrong generation.
  alphaSvalue     = settingsPtr->parm("SpaceShower:alphaSvalue");
  alphaSorder     = settingsPtr->mode("SpaceShower:alphaSorder");
  alphaS2pi       = 0.5 * alphaSvalue / M_PI;
  
  // Initialize alpha_strong generation.
  alphaS.init( alphaSvalue, alphaSorder); 
  
  // Lambda for 5, 4 and 3 flavours.
  Lambda5flav     = alphaS.Lambda5(); 
  Lambda4flav     = alphaS.Lambda4(); 
  Lambda3flav     = alphaS.Lambda3(); 
  Lambda5flav2    = pow2(Lambda5flav);
  Lambda4flav2    = pow2(Lambda4flav);
  Lambda3flav2    = pow2(Lambda3flav);
 
  // Regularization of QCD evolution for pT -> 0. Can be taken 
  // same as for multiparton interactions, or be set separately.
  useSamePTasMPI  = settingsPtr->flag("SpaceShower:samePTasMPI"); 
  if (useSamePTasMPI) {
    pT0Ref        = settingsPtr->parm("MultipartonInteractions:pT0Ref");
    ecmRef        = settingsPtr->parm("MultipartonInteractions:ecmRef");
    ecmPow        = settingsPtr->parm("MultipartonInteractions:ecmPow");
    pTmin         = settingsPtr->parm("MultipartonInteractions:pTmin");
  } else {
    pT0Ref        = settingsPtr->parm("SpaceShower:pT0Ref");
    ecmRef        = settingsPtr->parm("SpaceShower:ecmRef");
    ecmPow        = settingsPtr->parm("SpaceShower:ecmPow");
    pTmin         = settingsPtr->parm("SpaceShower:pTmin");
  }

  // Calculate nominal invariant mass of events. Set current pT0 scale.
  sCM             = m2( beamAPtr->p(), beamBPtr->p());
  eCM             = sqrt(sCM);
  pT0             = pT0Ref * pow(eCM / ecmRef, ecmPow);

  // Restrict pTmin to ensure that alpha_s(pTmin^2 + pT_0^2) does not blow up.
  double pTminAbs = sqrtpos(pow2(LAMBDA3MARGIN) * Lambda3flav2 / renormMultFac
                  - pT0*pT0);
  if (pTmin < pTminAbs) { 
    pTmin         = pTminAbs;
    ostringstream newPTmin;
    newPTmin << fixed << setprecision(3) << pTmin;
    infoPtr->errorMsg("Warning in SpaceShower::init: pTmin too low",
		      ", raised to " + newPTmin.str() );
    infoPtr->setTooLowPTmin(true);
  }

  // Parameters of alphaEM generation.
  alphaEMorder    = settingsPtr->mode("SpaceShower:alphaEMorder");

  // Initialize alphaEM generation.
  alphaEM.init( alphaEMorder, settingsPtr); 
 
  // Parameters of QED evolution.
  pTminChgQ       = settingsPtr->parm("SpaceShower:pTminchgQ"); 
  pTminChgL       = settingsPtr->parm("SpaceShower:pTminchgL"); 

  // Derived parameters of QCD evolution.
  pT20            = pow2(pT0);
  pT2min          = pow2(pTmin);
  pT2minChgQ      = pow2(pTminChgQ);
  pT2minChgL      = pow2(pTminChgL);

  // Various other parameters. 
  doMEcorrections = settingsPtr->flag("SpaceShower:MEcorrections");
  doMEafterFirst  = settingsPtr->flag("SpaceShower:MEafterFirst");
  doPhiPolAsym    = settingsPtr->flag("SpaceShower:phiPolAsym");
  doPhiIntAsym    = settingsPtr->flag("SpaceShower:phiIntAsym");
  strengthIntAsym = settingsPtr->parm("SpaceShower:strengthIntAsym");
  nQuarkIn        = settingsPtr->mode("SpaceShower:nQuarkIn");

  // Possibility of two predetermined hard emissions in event.
  doSecondHard    = settingsPtr->flag("SecondHard:generate");

  // Optional dampening at small pT's when large multiplicities.
  enhanceScreening 
    = settingsPtr->mode("MultipartonInteractions:enhanceScreening");
  if (!useSamePTasMPI) enhanceScreening = 0;

  // Possibility to allow user veto of emission step.
  canVetoEmission = (userHooksPtr != 0) 
                  ? userHooksPtr->canVetoISREmission() : false;

} 

//--------------------------------------------------------------------------

// Find whether to limit maximum scale of emissions.
// Also allow for dampening at factorization or renormalization scale. 

bool SpaceShower::limitPTmax( Event& event, double Q2Fac, double Q2Ren) {

  // Find whether to limit pT. Begin by user-set cases.
  bool dopTlimit = false;
  dopTlimit1 = dopTlimit2 = false;
  if      (pTmaxMatch == 1) dopTlimit = dopTlimit1 = dopTlimit2 = true;
  else if (pTmaxMatch == 2) dopTlimit = dopTlimit1 = dopTlimit2 = false;
   
  // Look if any quark (u, d, s, c, b), gluon or photon in final state. 
  else {
    int n21 = 0;
    for (int i = 5; i < event.size(); ++i) {
      if (event[i].status() == -21) ++n21;
      else if (n21 == 0) {
        int idAbs = event[i].idAbs();
        if (idAbs <= 5 || idAbs == 21 || idAbs == 22) dopTlimit1 = true;
      } else if (n21 == 2) {
        int idAbs = event[i].idAbs();
        if (idAbs <= 5 || idAbs == 21 || idAbs == 22) dopTlimit2 = true;
      }
    }
    dopTlimit = (doSecondHard) ? (dopTlimit1 && dopTlimit2) : dopTlimit1;  
  }

  // Dampening at factorization or renormalization scale; only for hardest.
  dopTdamp   = false;
  pT2damp    = 0.;
  if ( !dopTlimit1 && (pTdampMatch == 1 || pTdampMatch == 2) ) {
    dopTdamp = true;
    pT2damp  = pow2(pTdampFudge) * ((pTdampMatch == 1) ? Q2Fac : Q2Ren);
  }

  // Done.
  return dopTlimit;
 
}

//--------------------------------------------------------------------------

// Prepare system for evolution; identify ME.
// Routine may be called after multiparton interactions, for a new subystem.

void SpaceShower::prepare( int iSys, Event& event, bool limitPTmaxIn) {

  // Find positions of incoming colliding partons.
  int in1 = partonSystemsPtr->getInA(iSys);
  int in2 = partonSystemsPtr->getInB(iSys);

  // Rescattered partons cannot radiate.
  bool canRadiate1 = !(event[in1].isRescatteredIncoming());
  bool canRadiate2 = !(event[in2].isRescatteredIncoming());

  // Reset dipole-ends list for first interaction. Also resonances.
  if (iSys == 0) dipEnd.resize(0);
  if (iSys == 0) idResFirst  = 0;
  if (iSys == 1) idResSecond = 0;

  // Find matrix element corrections for system.
  int MEtype = findMEtype( iSys, event); 

  // In case of DPS overwrite limitPTmaxIn by saved value.
  if (doSecondHard && iSys == 0) limitPTmaxIn = dopTlimit1; 
  if (doSecondHard && iSys == 1) limitPTmaxIn = dopTlimit2; 

  // Maximum pT scale for dipole ends.
  double pTmax1 = (limitPTmaxIn) ? event[in1].scale() : eCM;
  double pTmax2 = (limitPTmaxIn) ? event[in2].scale() : eCM;
  if ( limitPTmaxIn && (iSys == 0 || (iSys == 1 && doSecondHard)) ) {
    pTmax1 *= pTmaxFudge;
    pTmax2 *= pTmaxFudge;
  } else if (limitPTmaxIn && iSys > 0) {
    pTmax1 *= pTmaxFudgeMPI;
    pTmax2 *= pTmaxFudgeMPI;
  }

  // Find dipole ends for QCD radiation.
  // Note: colour type can change during evolution, so book also if zero.
  if (doQCDshower) {
    int colType1 = event[in1].colType();
    if (canRadiate1) dipEnd.push_back( SpaceDipoleEnd( iSys,  1, 
      in1, in2, pTmax1, colType1, 0, MEtype, canRadiate2) );
    int colType2 = event[in2].colType();
    if (canRadiate2) dipEnd.push_back( SpaceDipoleEnd( iSys,  2, 
      in2, in1, pTmax2, colType2, 0, MEtype, canRadiate1) );
  }

  // Find dipole ends for QED radiation.
  // Note: charge type can change during evolution, so book also if zero.
  if (doQEDshowerByQ || doQEDshowerByL) {
    int chgType1 = ( (event[in1].isQuark() && doQEDshowerByQ)
      || (event[in1].isLepton() && doQEDshowerByL) )
      ? event[in1].chargeType() : 0;
    // Special: photons have charge zero, but can evolve (only off Q for now)
    if (event[in1].id() == 22 && doQEDshowerByQ) chgType1 = 22 ;
    if (canRadiate1) dipEnd.push_back( SpaceDipoleEnd( iSys, -1, 
      in1, in2, pTmax1, 0, chgType1, MEtype, canRadiate2) );
    int chgType2 = ( (event[in2].isQuark() && doQEDshowerByQ)
      || (event[in2].isLepton() && doQEDshowerByL) )
      ? event[in2].chargeType() : 0;
    // Special: photons have charge zero, but can evolve (only off Q for now)
    if (event[in2].id() == 22 && doQEDshowerByQ) chgType2 = 22 ;
    if (canRadiate2) dipEnd.push_back( SpaceDipoleEnd( iSys, -2, 
      in2, in1, pTmax2, 0, chgType2, MEtype, canRadiate1) );
  }

  // Store the z and pT2 values of the last previous splitting
  // when an event history has already been constructed.
  if (iSys == 0 && infoPtr->hasHistory()) {
    double zNow   = infoPtr->zNowISR();
    double pT2Now = infoPtr->pT2NowISR();
    for (int iDipEnd = 0; iDipEnd < int(dipEnd.size()); ++iDipEnd) {
      dipEnd[iDipEnd].zOld = zNow;
      dipEnd[iDipEnd].pT2Old = pT2Now;
      ++dipEnd[iDipEnd].nBranch;
    }
  }

}

//--------------------------------------------------------------------------
 
// Select next pT in downwards evolution of the existing dipoles.

double SpaceShower::pTnext( Event& event, double pTbegAll, double pTendAll, 
  int nRadIn) {

  // Current cm energy, in case it varies between events.
  sCM           = m2( beamAPtr->p(), beamBPtr->p());
  eCM           = sqrt(sCM);
  pTbegRef      = pTbegAll;

  // Starting values: no radiating dipole found.
  nRad          = nRadIn;
  double pT2sel = pow2(pTendAll);
  iDipSel       = 0;
  iSysSel       = 0;
  dipEndSel     = 0; 

  // Loop over all possible dipole ends.
  for (int iDipEnd = 0; iDipEnd < int(dipEnd.size()); ++iDipEnd) {
    iDipNow        = iDipEnd;
    dipEndNow      = &dipEnd[iDipEnd];        
    iSysNow        = dipEndNow->system;
    dipEndNow->pT2 = 0.;
   
    // Check whether dipole end should be allowed to shower. 
    double pT2begDip = pow2( min( pTbegAll, dipEndNow->pTmax ));
    if (pT2begDip > pT2sel 
      && ( dipEndNow->colType != 0 || dipEndNow->chgType != 0 ) ) {
      double pT2endDip = 0.;

      // Determine lower cut for evolution, for QCD or QED (q or l).      
      if (dipEndNow->colType != 0) pT2endDip = max( pT2sel, pT2min );   
      else if (abs(dipEndNow->chgType) != 3) pT2endDip 
        = max( pT2sel, pT2minChgQ );   
      else pT2endDip = max( pT2sel, pT2minChgL );  

      // Find properties of dipole and radiating dipole end.
      sideA        = ( abs(dipEndNow->side) == 1 ); 
      BeamParticle& beamNow = (sideA) ? *beamAPtr : *beamBPtr;
      BeamParticle& beamRec = (sideA) ? *beamBPtr : *beamAPtr;
      iNow         = beamNow[iSysNow].iPos();
      iRec         = beamRec[iSysNow].iPos();
      idDaughter   = beamNow[iSysNow].id();
      xDaughter    = beamNow[iSysNow].x();
      x1Now        = (sideA) ? xDaughter : beamRec[iSysNow].x();
      x2Now        = (sideA) ? beamRec[iSysNow].x() : xDaughter;

      // Note dipole mass correction when recoiler is a rescatter.
      m2Rec        = (dipEndNow->normalRecoil) ? 0. : event[iRec].m2(); 
      m2Dip        = x1Now * x2Now * sCM + m2Rec;

      // Now do evolution in pT2, for QCD or QED 
      if (pT2begDip > pT2endDip) { 
        if (dipEndNow->colType != 0) pT2nextQCD( pT2begDip, pT2endDip);
        else                         pT2nextQED( pT2begDip, pT2endDip);
      }

      // Update if found larger pT than current maximum.
      if (dipEndNow->pT2 > pT2sel) {
        pT2sel    = dipEndNow->pT2;
        iDipSel   = iDipNow;
        iSysSel   = iSysNow;
        dipEndSel = dipEndNow;
      }

    // End loop over dipole ends.
    }
  } 

  // Return nonvanishing value if found pT is bigger than already found.
  return (dipEndSel == 0) ? 0. : sqrt(pT2sel); 
}

//--------------------------------------------------------------------------

// Evolve a QCD dipole end. 

void SpaceShower::pT2nextQCD( double pT2begDip, double pT2endDip) { 

  // Some properties and kinematical starting values.
  BeamParticle& beam = (sideA) ? *beamAPtr : *beamBPtr;
  bool   isGluon     = (idDaughter == 21);
  bool   isValence   = beam[iSysNow].isValence();
  int    MEtype      = dipEndNow->MEtype;
  double pT2         = pT2begDip;
  double xMaxAbs     = beam.xMax(iSysNow);
  double zMinAbs     = xDaughter / xMaxAbs;
  if (xMaxAbs < 0.) {
    infoPtr->errorMsg("Warning in SpaceShower::pT2nextQCD: "
    "xMaxAbs negative"); 
    return;
  }

  // Starting values for handling of massive quarks (c/b), if any.
  double idMassive   = 0;
  if ( abs(idDaughter) == 4 ) idMassive = 4;
  if ( abs(idDaughter) == 5 ) idMassive = 5;
  bool   isMassive   = (idMassive > 0);
  double m2Massive   = 0.;
  double mRatio      = 0.;
  double zMaxMassive = 1.;
  double m2Threshold = pT2;

  // Evolution below scale of massive quark or at large x is impossible.
  if (isMassive) { 
    m2Massive = (idMassive == 4) ? m2c : m2b;
    if (pT2 < HEAVYPT2EVOL * m2Massive) return;
    mRatio = sqrt( m2Massive / m2Dip );
    zMaxMassive = (1. -  mRatio) / ( 1. +  mRatio * (1. -  mRatio) ); 
    if (xDaughter > HEAVYXEVOL * zMaxMassive * xMaxAbs) return; 
  
    // Find threshold scale below which only g -> Q + Qbar will be allowed.
    m2Threshold = (idMassive == 4) ? min( pT2, CTHRESHOLD * m2c)
      : min( pT2, BTHRESHOLD * m2b); 
  }
  
  // Variables used inside evolution loop. (Mainly dummy starting values.)
  int    nFlavour       = 3; 
  double b0             = 4.5;
  double Lambda2        = Lambda3flav2;
  double pT2minNow      = pT2endDip; 
  int    idMother       = 0; 
  int    idSister       = 0;
  double z              = 0.;
  double zMaxAbs        = 0.;
  double zRootMax       = 0.;
  double zRootMin       = 0.;
  double g2gInt         = 0.; 
  double q2gInt         = 0.; 
  double q2qInt         = 0.;
  double g2qInt         = 0.;
  double g2Qenhance     = 0.;
  double xPDFdaughter   = 0.;
  double xPDFmother[21] = {0.};
  double xPDFgMother    = 0.;
  double xPDFmotherSum  = 0.;
  double kernelPDF      = 0.;
  double xMother        = 0.;
  double wt             = 0.;
  double Q2             = 0.;
  double mSister        = 0.;
  double m2Sister       = 0.;
  double pT2corr        = 0.;
  double pT2PDF         = pT2;
  bool   needNewPDF     = true;

  // Begin evolution loop towards smaller pT values.
  int    loopTinyPDFdau = 0;
  bool   hasTinyPDFdau  = false;
  do { 
    wt = 0.;

    // Bad sign if repeated looping with small daughter PDF, so fail.
    // (Example: if all PDF's = 0 below Q_0, except for c/b companion.)
    if (hasTinyPDFdau) ++loopTinyPDFdau;  
    if (loopTinyPDFdau > MAXLOOPTINYPDF) {
      infoPtr->errorMsg("Warning in SpaceShower::pT2nextQCD: "
      "small daughter PDF"); 
      return;
    }

    // Initialize integrals of splitting kernels and evaluate parton 
    // densities at the beginning. Reinitialize after long evolution 
    // in pT2 or when crossing c and b flavour thresholds.
    if (needNewPDF || pT2 < EVALPDFSTEP * pT2PDF) {
      pT2PDF        = pT2;
      hasTinyPDFdau = false;

      // Determine overestimated z range; switch at c and b masses.
      if (pT2 > m2b) {
        nFlavour  = 5;
        pT2minNow = max( m2b, pT2endDip);
        b0        = 23./6.;
        Lambda2   = Lambda5flav2;
      } else if (pT2 > m2c) {
        nFlavour  = 4;
        pT2minNow = max( m2c, pT2endDip);
        b0        = 25./6.;
        Lambda2   = Lambda4flav2;
      } else { 
        nFlavour  = 3;
        pT2minNow = pT2endDip;
        b0        = 27./6.;
        Lambda2   = Lambda3flav2;
      }
      // A change of renormalization scale expressed by a change of Lambda. 
      Lambda2    /= renormMultFac;
      zMaxAbs     = 1. - 0.5 * (pT2minNow / m2Dip) *
        ( sqrt( 1. + 4. * m2Dip / pT2minNow ) - 1. );
      if (isMassive) zMaxAbs = min( zMaxAbs, zMaxMassive); 

      // Go to another z range with lower mass scale if current is closed.
      if (zMinAbs > zMaxAbs) { 
        if (nFlavour == 3 || (idMassive == 4 && nFlavour == 4) 
          || idMassive == 5) return;
        pT2 = (nFlavour == 4) ? m2c : m2b;
        continue;
      } 

      // Parton density of daughter at current scale. 
      xPDFdaughter = beam.xfISR(iSysNow, idDaughter, xDaughter, 
        factorMultFac * pT2);
      if (xPDFdaughter < TINYPDF) {
        xPDFdaughter  = TINYPDF;
        hasTinyPDFdau = true;
      }

      // Integrals of splitting kernels for gluons: g -> g, q -> g.
      if (isGluon) {
        g2gInt = 6. * log(zMaxAbs * (1.-zMinAbs) 
          / (zMinAbs * (1.-zMaxAbs)));
        if (doMEcorrections) g2gInt *= calcMEmax(MEtype, 21, 21);
        q2gInt = (16./3.) * (1./sqrt(zMinAbs) - 1./sqrt(zMaxAbs));
        if (doMEcorrections) q2gInt *= calcMEmax(MEtype, 1, 21);

        // Parton density of potential quark mothers to a g.
        xPDFmotherSum = 0.;
        for (int i = -nQuarkIn; i <= nQuarkIn; ++i) {
          if (i == 0) {
            xPDFmother[10] = 0.;
          } else {
            xPDFmother[i+10] = beam.xfISR(iSysNow, i, xDaughter, 
              factorMultFac * pT2); 
            xPDFmotherSum += xPDFmother[i+10]; 
          }
        } 

        // Total QCD evolution coefficient for a gluon.
        kernelPDF = g2gInt + q2gInt * xPDFmotherSum / xPDFdaughter;

      // For valence quark only need consider q -> q g branchings.
      // Introduce an extra factor sqrt(z) to smooth bumps.
      } else if (isValence) {
        zRootMin = (1. + sqrt(zMinAbs)) / (1. - sqrt(zMinAbs));
        zRootMax = (1. + sqrt(zMaxAbs)) / (1. - sqrt(zMaxAbs));
        q2qInt = (8./3.) * log( zRootMax / zRootMin );
        if (doMEcorrections) q2qInt *= calcMEmax(MEtype, 1, 1);
        kernelPDF = q2qInt; 

      // Integrals of splitting kernels for quarks: q -> q, g -> q.
      } else {
        q2qInt = (8./3.) * log( (1. - zMinAbs) / (1. - zMaxAbs) );
        if (doMEcorrections) q2qInt *= calcMEmax(MEtype, 1, 1);
        g2qInt = 0.5 * (zMaxAbs - zMinAbs);
        if (doMEcorrections) g2qInt *= calcMEmax(MEtype, 21, 1);

        // Increase estimated upper weight for g -> Q + Qbar.
        if (isMassive) {
          if (alphaSorder == 0) g2Qenhance = log(pT2/m2Massive) 
            / log(m2Threshold/m2Massive);    
          else {
            double m2log = log( m2Massive / Lambda2);
            g2Qenhance = log( log(pT2/Lambda2) / m2log ) 
              / log( log(m2Threshold/Lambda2) / m2log );
          }
          g2qInt *= g2Qenhance;
        }

        // Parton density of a potential gluon mother to a q.
        xPDFgMother = beam.xfISR(iSysNow, 21, xDaughter, factorMultFac * pT2);

        // Total QCD evolution coefficient for a quark.
        kernelPDF = q2qInt + g2qInt * xPDFgMother / xPDFdaughter;
      }

      // End evaluation of splitting kernels and parton densities.
      needNewPDF = false;
    }
    if (kernelPDF < TINYKERNELPDF) return;

    // Pick pT2 (in overestimated z range), for one of three different cases.
    // Assume form alphas(pT0^2 + pT^2) * dpT^2/(pT0^2 + pT^2).
    double Q2alphaS;

    // Fixed alpha_strong.
    if (alphaSorder == 0) {
      pT2 = (pT2 + pT20) * pow( rndmPtr->flat(), 
        1. / (alphaS2pi * kernelPDF)) - pT20;

    // First-order alpha_strong.
    } else if (alphaSorder == 1) {
      pT2 = Lambda2 * pow( (pT2 + pT20) / Lambda2, 
        pow(rndmPtr->flat(), b0 / kernelPDF) ) - pT20;

    // For second order reject by second term in alpha_strong expression.
    } else {
      do {
        pT2 = Lambda2 * pow( (pT2 + pT20) / Lambda2,
          pow(rndmPtr->flat(), b0 / kernelPDF) ) - pT20;
        Q2alphaS = renormMultFac * max( pT2 + pT20, 
          pow2(LAMBDA3MARGIN) * Lambda3flav2);
      } while (alphaS.alphaS2OrdCorr(Q2alphaS) < rndmPtr->flat()
        && pT2 > pT2minNow);
    }

    // Check for pT2 values that prompt special action.

    // If fallen into b threshold region, force g -> b + bbar.
    if (idMassive == 5 && pT2 < m2Threshold) {
      pT2nearQCDthreshold( beam, m2Massive, m2Threshold, xMaxAbs, 
        zMinAbs, zMaxMassive );
      return;

    // If crossed b threshold, continue evolution from this threshold.
    } else if (nFlavour == 5 && pT2 < m2b) {  
      needNewPDF = true;
      pT2 = m2b;
      continue;

    // If fallen into c threshold region, force g -> c + cbar.
    } else if (idMassive == 4 && pT2 < m2Threshold) {
      pT2nearQCDthreshold( beam, m2Massive, m2Threshold, xMaxAbs, 
        zMinAbs, zMaxMassive );
      return; 

    // If crossed c threshold, continue evolution from this threshold.
    } else if (nFlavour == 4 && pT2 < m2c) { 
      needNewPDF = true;
      pT2 = m2c;
      continue;

    // Abort evolution if below cutoff scale, or below another branching.
    } else if (pT2 < pT2endDip) return; 

    // Select z value of branching to g, and corrective weight.
    if (isGluon) {
      // g -> g (+ g). 
      if (rndmPtr->flat() * kernelPDF < g2gInt) {
        idMother = 21;
        idSister = 21;
        z = 1. / ( 1. + ((1. - zMinAbs) / zMinAbs) * pow( (zMinAbs * 
          (1. - zMaxAbs)) / (zMaxAbs * (1. - zMinAbs)), rndmPtr->flat() ) );
        wt = pow2( 1. - z * (1. - z));
      } else {
      // q -> g (+ q): also select flavour. 
        double temp = xPDFmotherSum * rndmPtr->flat();
        idMother = -nQuarkIn - 1;
        do { temp -= xPDFmother[(++idMother) + 10]; } 
        while (temp > 0. && idMother < nQuarkIn);  
        idSister = idMother;
        z = (zMinAbs * zMaxAbs) / pow2( sqrt(zMinAbs) + rndmPtr->flat() 
          * ( sqrt(zMaxAbs)- sqrt(zMinAbs) ));
        wt = 0.5 * (1. + pow2(1. - z)) * sqrt(z) 
          * xPDFdaughter / xPDFmother[idMother + 10];
      } 

    // Select z value of branching to q, and corrective weight.
    // Include massive kernel corrections for c and b quarks.
    } else {
      // q -> q (+ g). 
      if (isValence || rndmPtr->flat() * kernelPDF < q2qInt) {
        idMother = idDaughter;
        idSister = 21;
        // Valence more peaked at large z.
        if (isValence) {
          double zTmp = zRootMin * pow(zRootMax / zRootMin, rndmPtr->flat() );
          z = pow2( (1. - zTmp) / (1. + zTmp) );
        } else {
          z = 1. - (1. - zMinAbs) * pow( (1. - zMaxAbs) / (1. - zMinAbs),
            rndmPtr->flat() );
        } 
        if (!isMassive) { 
          wt = 0.5 * (1. + pow2(z));
        } else {
        //?? Bug?? should be 2 more for massive part??
        //  wt = 0.5 * (1. + pow2(z) - z * pow2(1.-z) * m2Massive / pT2);
          wt = 0.5 * (1. + pow2(z)) - z * pow2(1.-z) * m2Massive / pT2;
        }
        if (isValence) wt *= sqrt(z);
      // g -> q (+ qbar). 
      } else {
        idMother = 21;
        idSister = - idDaughter; 
        z = zMinAbs + rndmPtr->flat() * (zMaxAbs - zMinAbs);
        if (!isMassive) { 
          wt = (pow2(z) + pow2(1.-z)) * xPDFdaughter / xPDFgMother ;
        } else {
          wt = (pow2(z) + pow2(1.-z) + 2. * z * (1.-z) * m2Massive / pT2) 
            * xPDFdaughter / (xPDFgMother * g2Qenhance) ;
        }
      }
    }

    // Derive Q2 and x of mother from pT2 and z. 
    Q2      = pT2 / (1. - z);
    xMother = xDaughter / z;
    // Correction to x for massive recoiler from rescattering.
    if (!dipEndNow->normalRecoil) {
      if (sideA) xMother += (m2Rec / (x2Now * sCM)) * (1. / z - 1.);
      else       xMother += (m2Rec / (x1Now * sCM)) * (1. / z - 1.);
    }
    if(xMother > xMaxAbs) { wt = 0.; continue; }

    // Forbidden emission if outside allowed z range for given pT2.
    mSister = particleDataPtr->m0(idSister);
    m2Sister = pow2(mSister);
    pT2corr = Q2 - z * (m2Dip + Q2) * (Q2 + m2Sister) / m2Dip;
    if(pT2corr < TINYPT2) { wt = 0.; continue; }

    // Optionally veto emissions not ordered in rapidity (= angle).
    if ( doRapidityOrder && dipEndNow->nBranch > 0
      && pT2 > pow2( (1. - z) / (z * (1. - dipEndNow->zOld)) ) 
      * dipEndNow->pT2Old ) { wt = 0.; continue; }

    // If creating heavy quark by Q -> g + Q then next need g -> Q + Qbar.
    // So minimum total mass2 is 4 * m2Sister, but use more to be safe.
    if ( isGluon && ( abs(idMother) == 4 || abs(idMother) == 5 )) {
      double m2QQsister =  EXTRASPACEQ * 4. * m2Sister;
      double pT2QQcorr = Q2 - z * (m2Dip + Q2) * (Q2 + m2QQsister) / m2Dip;
      if(pT2QQcorr < TINYPT2) { wt = 0.; continue; }
    }  

    // Evaluation of ME correction.
    if (doMEcorrections) wt *= calcMEcorr(MEtype, idMother, idDaughter, 
      m2Dip, z, Q2) / calcMEmax(MEtype, idMother, idDaughter); 

    // Optional dampening of large pT values in first radiation.
    if (dopTdamp && iSysNow == 0 && MEtype == 0 && nRad == 0) 
      wt *= pT2damp / (pT2 + pT2damp);

    // Idea suggested by Gosta Gustafson: increased screening in events
    // with large activity can be simulated by pT0_eff = sqrt(n) * pT0. 
    if (enhanceScreening == 2) {
      int nSysNow     = infoPtr->nMPI() + infoPtr->nISR() + 1;
      double WTscreen = pow2( (pT2 + pT20) / (pT2 + nSysNow * pT20) );
      wt             *= WTscreen;
    } 

    // Evaluation of new daughter and mother PDF's.
    double xPDFdaughterNew = max ( TINYPDF, 
      beam.xfISR(iSysNow, idDaughter, xDaughter, factorMultFac * pT2) );
    double xPDFmotherNew = 
      beam.xfISR(iSysNow, idMother, xMother, factorMultFac * pT2);
    wt *= xPDFmotherNew / xPDFdaughterNew;

    // Check that valence step does not cause problem.
    if (wt > 1. && pT2 > PT2MINWARN) infoPtr->errorMsg("Warning in "
      "SpaceShower::pT2nextQCD: weight above unity"); 

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < rndmPtr->flat()) ;

  // Save values for (so far) acceptable branching.
  dipEndNow->store( idDaughter,idMother, idSister, x1Now, x2Now, m2Dip,
    pT2, z, xMother, Q2, mSister, m2Sister, pT2corr);  

}

//--------------------------------------------------------------------------

// Evolve a QCD dipole end near threshold, with g -> Q + Qbar enforced.
// Note: No explicit Sudakov factor formalism here. Instead use that 
// df_Q(x, pT2) = (alpha_s/2pi) * (dT2/pT2) * ((gluon) * (splitting)).
// This implies that effects of Q -> Q + g are neglected in this range. 

void SpaceShower::pT2nearQCDthreshold( BeamParticle& beam, 
  double m2Massive, double m2Threshold, double xMaxAbs, 
  double zMinAbs, double zMaxMassive) {

  // Initial values, to be used in kinematics and weighting.
  double Lambda2       = (abs(idDaughter) == 4) ? Lambda4flav2 : Lambda5flav2;
  Lambda2             /= renormMultFac;
  double logM2Lambda2  = log( m2Massive / Lambda2 );
  double xPDFmotherOld = beam.xfISR(iSysNow, 21, xDaughter, 
    factorMultFac * m2Threshold);
  if (xPDFmotherOld == 0) {
    infoPtr->errorMsg("Error in SpaceShower::pT2nearQCDthreshold: "
      "xPDFmotherOld is 0"); 
    return; 
  }

  // Variables used inside evolution loop. (Mainly dummy start values.)
  int    loop    = 0;
  double wt      = 0.;
  double pT2     = 0.; 
  double z       = 0.; 
  double Q2      = 0.; 
  double pT2corr = 0.;
  double xMother = 0.;

  // Begin loop over tries to find acceptable g -> Q + Qbar branching. 
  do { 
    wt = 0.;

    // Check that not caught in infinite loop with impossible kinematics.
    if (++loop > 100) { 
      infoPtr->errorMsg("Error in SpaceShower::pT2nearQCDthreshold: "
        "stuck in loop"); 
      return; 
    }

    // Pick dpT2/pT2 in range [m2Massive,thresholdRatio * m2Massive]. 
    pT2 = m2Massive * pow( m2Threshold / m2Massive, rndmPtr->flat() ); 

    // Pick z flat in allowed range.
    z = zMinAbs + rndmPtr->flat() * (zMaxMassive - zMinAbs);

    // Check that kinematically possible choice.
    Q2 = pT2 / (1.-z) - m2Massive;
    pT2corr = Q2 - z * (m2Dip + Q2) * (Q2 + m2Massive) / m2Dip;
    if(pT2corr < TINYPT2) continue;
    
    // Correction factor for running alpha_s.  ??
    wt = logM2Lambda2 / log( pT2 / Lambda2 ); 

    // Correction factor for splitting kernel.
    wt *= pow2(z) + pow2(1.-z) + 2. * z * (1.-z) * m2Massive / pT2;

    // x, including correction for massive recoiler from rescattering.
    xMother = xDaughter / z;
    if (!dipEndNow->normalRecoil) {
      if (sideA) xMother += (m2Rec / (x2Now * sCM)) * (1. / z - 1.);
      else       xMother += (m2Rec / (x1Now * sCM)) * (1. / z - 1.);
    }
    if (xMother > xMaxAbs) { wt = 0.; continue; }

    // Correction factor for gluon density.
    double xPDFmotherNew = beam.xfISR(iSysNow, 21, xMother, 
      factorMultFac * pT2);
    wt *= xPDFmotherNew / xPDFmotherOld;

  // Iterate until acceptable pT and z.
  } while (wt < rndmPtr->flat()) ;

  // Save values for (so far) acceptable branching.
  double mSister = (abs(idDaughter) == 4) ? mc : mb;  
  dipEndNow->store( idDaughter, 21, -idDaughter, x1Now, x2Now, m2Dip,
    pT2, z, xMother, Q2, mSister, pow2(mSister), pT2corr);  

}

//--------------------------------------------------------------------------

// Evolve a QED dipole end. 

void SpaceShower::pT2nextQED( double pT2begDip, double pT2endDip) { 

  // Type of dipole and starting values.
  BeamParticle& beam  = (sideA) ? *beamAPtr : *beamBPtr;
  bool   isLeptonBeam = beam.isLepton();
  int    MEtype       = dipEndNow->MEtype;
  bool   isPhoton     = (idDaughter == 22);
  double pT2          = pT2begDip;
  double m2Lepton     = (isLeptonBeam) ? pow2(beam.m()) : 0.; 
  if (isLeptonBeam && pT2begDip < m2Lepton) return;

  // Currently no f -> gamma branching implemented for lepton beams
  if (isPhoton && isLeptonBeam) return;

  // alpha_em at maximum scale provides upper estimate.
  double alphaEMmax  = alphaEM.alphaEM(renormMultFac * pT2begDip);
  double alphaEM2pi  = alphaEMmax / (2. * M_PI);

  // Maximum x of mother implies minimum z = xDaughter / xMother.
  double xMaxAbs  = (isLeptonBeam) ? LEPTONXMAX : beam.xMax(iSysNow);
  double zMinAbs  = xDaughter / xMaxAbs;
  if (xMaxAbs < 0.) {
    infoPtr->errorMsg("Warning in SpaceShower::pT2nextQED: "
    "xMaxAbs negative"); 
    return;
  }

  // Maximum z from minimum pT and, for lepton, from minimum x_gamma.
  double zMaxAbs = 1. - 0.5 * (pT2endDip / m2Dip) *
    ( sqrt( 1. + 4. * m2Dip / pT2endDip ) - 1. );
  if (isLeptonBeam) {
    double zMaxLepton = xDaughter / (xDaughter + LEPTONXMIN);
    if (zMaxLepton < zMaxAbs) zMaxAbs = zMaxLepton;
  }
  if (zMaxAbs < zMinAbs) return;

  // Variables used inside evolution loop. (Mainly dummy start values.)
  int    idMother = 0;
  int    idSister =22;
  double z        = 0.; 
  double xMother  = 0.; 
  double wt       = 0.; 
  double Q2       = 0.;
  double mSister  = 0.; 
  double m2Sister = 0.;
  double pT2corr  = 0.;
  
  // QED evolution of fermions
  if (!isPhoton) {

    // Integrals of splitting kernels for fermions: f -> f. Use 1 + z^2 < 2. 
    // Ansatz f(z) = 2 / (1 - z), with + 2 / (z - xDaughter) for lepton.
    double f2fInt  = 0.;
    double f2fIntA = 2. * log( (1. - zMinAbs) / (1. - zMaxAbs) );
    double f2fIntB = 0.;
    if (isLeptonBeam) {
      f2fIntB      = 2. * log( (zMaxAbs - xDaughter) / (zMinAbs - xDaughter) );
      f2fInt       = f2fIntA + f2fIntB; 
    } else f2fInt  = pow2(dipEndNow->chgType / 3.) * f2fIntA;
    
    // Upper estimate for evolution equation, including fudge factor. 
    if (doMEcorrections) f2fInt *= calcMEmax(MEtype, 1, 1);
    double kernelPDF = alphaEM2pi * f2fInt;
    double fudge = (isLeptonBeam) ? LEPTONFUDGE * log(m2Dip/m2Lepton) : 1.;
    kernelPDF *= fudge;
    if (kernelPDF < TINYKERNELPDF) return;
    
    // Begin evolution loop towards smaller pT values.
    do { 
      
      // Pick pT2 (in overestimated z range).
      // For l -> l gamma include extrafactor 1 / ln(pT2 / m2l) in evolution.
      double shift = pow(rndmPtr->flat(), 1. / kernelPDF);
      if (isLeptonBeam) pT2 = m2Lepton * pow( pT2 / m2Lepton, shift);
      else              pT2 = pT2 * shift; 
      
      // Abort evolution if below cutoff scale, or below another branching.
      if (pT2 < pT2endDip) return; 
      if (isLeptonBeam && pT2 < LEPTONPT2MIN * m2Lepton) return; 
      
      // Select z value of branching f -> f + gamma, and corrective weight.
      idMother = idDaughter;
      wt = 1.;
      if (isLeptonBeam) {
        if (f2fIntA > rndmPtr->flat() * (f2fIntA + f2fIntB)) { 
          z = 1. - (1. - zMinAbs) 
            * pow( (1. - zMaxAbs) / (1. - zMinAbs), rndmPtr->flat() );
        } else {
          z = xDaughter + (zMinAbs - xDaughter) 
            * pow( (zMaxAbs - xDaughter) / (zMinAbs - xDaughter), 
                  rndmPtr->flat() );  
        }
        wt *= (z - xDaughter) / (1. - xDaughter); 
      } else {
        z = 1. - (1. - zMinAbs) 
          * pow( (1. - zMaxAbs) / (1. - zMinAbs), rndmPtr->flat() ); 
      }
      wt *= 0.5 * (1. + pow2(z));
      
      // Derive Q2 and x of mother from pT2 and z. 
      Q2      = pT2 / (1. - z);
      xMother = xDaughter / z;
      // Correction to x for massive recoiler from rescattering.
      if (!dipEndNow->normalRecoil) {
        if (sideA) xMother += (m2Rec / (x2Now * sCM)) * (1. / z - 1.);
        else       xMother += (m2Rec / (x1Now * sCM)) * (1. / z - 1.);
      }
      if(xMother > xMaxAbs) { wt = 0.; continue; }
      
      // Forbidden emission if outside allowed z range for given pT2.
      mSister  = 0.;
      m2Sister = 0.;
      pT2corr  = Q2 - z * (m2Dip + Q2) * (Q2 + m2Sister) / m2Dip;
      if(pT2corr < TINYPT2) { wt = 0.; continue; }
      
      // Correct by ln(pT2 / m2l) and fudge factor.  
      if (isLeptonBeam) wt *= log(pT2 / m2Lepton) / fudge;
      
      // Evaluation of ME correction.
      if (doMEcorrections) wt *= calcMEcorr(MEtype, idMother, idDaughter, 
         m2Dip, z, Q2) / calcMEmax(MEtype, idMother, idDaughter);

      // Extra QED correction for f fbar -> W+- gamma. Debug??
      if (doMEcorrections && MEtype == 1 && idDaughter == idMother
        && ( (iSysNow == 0 && idResFirst  == 24)
          || (iSysNow == 1 && idResSecond == 24) ) ) {
        double tHe  = -Q2;
        double uHe  = Q2 - m2Dip * (1. - z) / z;
        double chg1 = abs(dipEndNow->chgType / 3.);
        double chg2 = 1. - chg1;
        wt *= pow2(chg1 * uHe - chg2 * tHe) 
          / ( (tHe + uHe) * (pow2(chg1) * uHe + pow2(chg2) * tHe) );  
      }

      // Optional dampening of large pT values in first radiation.
      if (dopTdamp && iSysNow == 0 && MEtype == 0 && nRad == 0) 
        wt *= pT2damp / (pT2 + pT2damp);

      // Correct to current value of alpha_EM.
      double alphaEMnow = alphaEM.alphaEM(renormMultFac * pT2);
      wt *= (alphaEMnow / alphaEMmax);
      
      // Evaluation of new daughter and mother PDF's.
      double xPDFdaughterNew = max ( TINYPDF, 
         beam.xfISR(iSysNow, idDaughter, xDaughter, factorMultFac * pT2) );
      double xPDFmotherNew   = 
         beam.xfISR(iSysNow, idMother, xMother, factorMultFac * pT2);
      wt *= xPDFmotherNew / xPDFdaughterNew;

    // Iterate until acceptable pT (or have fallen below pTmin).
    } while (wt < rndmPtr->flat()) ;
  }

  // QED evolution of photons (so far only for hadron beams).
  else {
    
    // Initial values
    int    nFlavour       = 3;         
    double kernelPDF      = 0.0;
    double xPDFdaughter   = 0.;
    double xPDFmother[21] = {0.};
    double xPDFmotherSum  = 0.0;
    double pT2PDF         = pT2;
    double pT2minNow      = pT2endDip;
    bool   needNewPDF     = true;

    // Begin evolution loop towards smaller pT values.
    int    loopTinyPDFdau = 0;
    bool   hasTinyPDFdau  = false;
    do { 
      wt = 0.;
      
      // Bad sign if repeated looping with small daughter PDF, so fail.
      if (hasTinyPDFdau) ++loopTinyPDFdau;  
      if (loopTinyPDFdau > MAXLOOPTINYPDF) {
        infoPtr->errorMsg("Warning in SpaceShower::pT2nextQED: "
                          "small daughter PDF"); 
        return;
      }

      // Initialize integrals of splitting kernels and evaluate parton 
      // densities at the beginning. Reinitialize after long evolution 
      // in pT2 or when crossing c and b flavour thresholds.
      if (needNewPDF || pT2 < EVALPDFSTEP * pT2PDF) {

        pT2PDF        = pT2;
        hasTinyPDFdau = false;

        // Determine overestimated z range; switch at c and b masses.
        if (pT2 > m2b && nQuarkIn >= 5) {
          nFlavour  = 5;
          pT2minNow = max( m2b, pT2endDip);
        } else if (pT2 > m2c && nQuarkIn >= 4) {
          nFlavour  = 4;
          pT2minNow = max( m2c, pT2endDip);
        } else { 
          nFlavour  = 3;
          pT2minNow = pT2endDip;
        }
        
        // Compute upper z limit 
        zMaxAbs = 1. - 0.5 * (pT2minNow / m2Dip) *
          ( sqrt( 1. + 4. * m2Dip / pT2minNow ) - 1. );

        // Parton density of daughter at current scale. 
        xPDFdaughter = beam.xfISR(iSysNow, idDaughter, xDaughter, 
          factorMultFac * pT2);
        if (xPDFdaughter < TINYPDF) {
          xPDFdaughter  = TINYPDF;
          hasTinyPDFdau = true;
        }
        
        // Integral over f -> gamma f splitting kernel.
        // Normalized so: 4/3 aS/2pi P(z) -> eq^2 * aEM/2pi P(z).
        // (Charge-weighting happens below.)
        double q2gInt = 4. * (1./sqrt(zMinAbs) - 1./sqrt(zMaxAbs));
        
        // Charge-weighted Parton density of potential quark mothers.
        xPDFmotherSum = 0.;
        for (int i = -nFlavour; i <= nFlavour; ++i) {
          if (i == 0) {
            xPDFmother[10] = 0.;
          } else {
            xPDFmother[i+10] = pow2((abs(i+1) % 2 + 1)/3.0) 
              * beam.xfISR(iSysNow, i, xDaughter, factorMultFac * pT2); 
            xPDFmotherSum += xPDFmother[i+10]; 
          }
        } 
        
        // Total QED evolution coefficient for a photon.
        kernelPDF = q2gInt * xPDFmotherSum / xPDFdaughter;
        
        // End evaluation of splitting kernels and parton densities.
        needNewPDF = false;
      }
      if (kernelPDF < TINYKERNELPDF) return;
      
      // Select pT2 for next trial branching 
      pT2 *= pow( rndmPtr->flat(), 1. / (alphaEM2pi * kernelPDF));

      // If crossed b threshold, continue evolution from this threshold.
      if (nFlavour == 5 && pT2 < m2b) {  
        needNewPDF = true;
        pT2 = m2b;
        continue;
      }

      // If crossed c threshold, continue evolution from this threshold.
      else if (nFlavour == 4 && pT2 < m2c) { 
        needNewPDF = true;
        pT2 = m2c;
        continue;
      }

      // Abort evolution if below cutoff scale, or below another branching.
      else if (pT2 < pT2endDip) return; 

      // Select flavour for trial branching
      double temp = xPDFmotherSum * rndmPtr->flat();
      idMother = -nQuarkIn - 1;
      do { 
        temp -= xPDFmother[(++idMother) + 10]; 
      } while (temp > 0. && idMother < nQuarkIn);  

      // Sister is same as mother, but can have m2 > 0
      idSister = idMother;
      mSister = particleDataPtr->m0(idSister);
      m2Sister = pow2(mSister);
      
      // Select z for trial branching
      z = (zMinAbs * zMaxAbs) / pow2( sqrt(zMinAbs) + rndmPtr->flat() 
                                      * ( sqrt(zMaxAbs)- sqrt(zMinAbs) ));
      
      // Trial weight: splitting kernel
      wt = 0.5 * (1. + pow2(1. - z)) * sqrt(z);
      
      // Trial weight: running alpha_EM
      double alphaEMnow = alphaEM.alphaEM(renormMultFac * pT2);
      wt *= (alphaEMnow / alphaEMmax);
      
      // Derive Q2 and x of mother from pT2 and z
      Q2      = pT2 / (1. - z);
      xMother = xDaughter / z;
      // Correction to x for massive recoiler from rescattering.
      if (!dipEndNow->normalRecoil) {
        if (sideA) xMother += (m2Rec / (x2Now * sCM)) * (1. / z - 1.);
        else       xMother += (m2Rec / (x1Now * sCM)) * (1. / z - 1.);
      }
      
      // Compute pT2corr
      pT2corr  = Q2 - z * (m2Dip + Q2) * (Q2 + m2Sister) / m2Dip;
      if(pT2corr < TINYPT2) { wt = 0.; continue; }
      
      // If creating heavy quark by Q -> gamma + Q then next g -> Q + Qbar.
      // So minimum total mass2 is 4 * m2Sister, but use more to be safe.
      if ( abs(idMother) == 4 || abs(idMother) == 5 ) {
        double m2QQsister =  EXTRASPACEQ * 4. * m2Sister;
        double pT2QQcorr = Q2 - z * (m2Dip + Q2) * (Q2 + m2QQsister) / m2Dip;
        if(pT2QQcorr < TINYPT2) { wt = 0.; continue; }
      }  
      
      // Optional dampening of large pT values in first radiation.
      if (dopTdamp && iSysNow == 0 && MEtype == 0 && nRad == 0) 
        wt *= pT2damp / (pT2 + pT2damp);

      // Evaluation of new daughter PDF 
      double xPDFdaughterNew = beam.xfISR(iSysNow, idDaughter, xDaughter, 
        factorMultFac * pT2);
      if (xPDFdaughterNew < TINYPDF) {
        xPDFdaughterNew = TINYPDF;
      }
      
      // Evaluation of new charge-weighted mother PDF 
      double xPDFmotherNew = pow2( (abs(idMother+1) % 2 + 1)/3.0 ) 
        * beam.xfISR(iSysNow, idMother, xMother, factorMultFac * pT2);
      
      // Trial weight: divide out old pdf ratio
      wt *= xPDFdaughter / xPDFmother[idMother + 10];
      
      // Trial weight: new pdf ratio
      wt *= xPDFmotherNew / xPDFdaughterNew;

    // Iterate until acceptable pT (or have fallen below pTmin).
    } while (wt < rndmPtr->flat()) ;    
  }

  // Save values for (so far) acceptable branching.
  dipEndNow->store( idDaughter, idMother, idSister, x1Now, x2Now, m2Dip,
    pT2, z, xMother, Q2, mSister, m2Sister, pT2corr);  

}

//--------------------------------------------------------------------------

// Kinematics of branching.
// Construct mother -> daughter + sister, with recoiler on other side. 

bool SpaceShower::branch( Event& event) {
  
  // Side on which branching occured.
  int side          = abs(dipEndSel->side);
  double sideSign   = (side == 1) ? 1. : -1.;

  // Read in flavour and colour variables.
  int iDaughter     = partonSystemsPtr->getInA(iSysSel);
  int iRecoiler     = partonSystemsPtr->getInB(iSysSel);
  if (side == 2) swap(iDaughter, iRecoiler);
  int idDaughterNow = dipEndSel->idDaughter;
  int idMother      = dipEndSel->idMother;
  int idSister      = dipEndSel->idSister;
  int colDaughter   = event[iDaughter].col();
  int acolDaughter  = event[iDaughter].acol();

  // Recoil parton may be rescatterer, requiring special processing. 
  bool normalRecoil = dipEndSel->normalRecoil; 
  int iRecoilMother = event[iRecoiler].mother1();

  // Read in kinematical variables.
  double x1         = dipEndSel->x1;
  double x2         = dipEndSel->x2;
  double xMo        = dipEndSel->xMo;
  double m2         = dipEndSel->m2Dip;
  double m          = sqrt(m2);
  double pT2        = dipEndSel->pT2;
  double z          = dipEndSel->z;
  double Q2         = dipEndSel->Q2; 
  double mSister    = dipEndSel->mSister;
  double m2Sister   = dipEndSel->m2Sister;
  double pT2corr    = dipEndSel->pT2corr;
  double x1New      = (side == 1) ? xMo : x1;
  double x2New      = (side == 2) ? xMo : x2;

  // Rescatter: kinematics may fail; use the rescatterFail flag to tell
  //            parton level to try again.
  rescatterFail     = false;

  // Construct kinematics of mother, sister and recoiler in old rest frame.
  // Normally both mother and recoiler are taken massless.
  double eNewRec, pzNewRec, pTbranch, pzMother;
  if (normalRecoil) {
    eNewRec         = 0.5 * (m2 + Q2) / m;
    pzNewRec        = -sideSign * eNewRec;
    pTbranch        = sqrt(pT2corr) * m2 / ( z * (m2 + Q2) );
    pzMother        = sideSign * 0.5 * m * ( (m2 - Q2) / ( z * (m2 + Q2) )
                    + (Q2 + m2Sister) / m2 ); 
  // More complicated kinematics when recoiler not massless. May fail.
  } else {
    m2Rec           = event[iRecoiler].m2();
    double s1Tmp    = m2 + Q2 - m2Rec;
    double s3Tmp    = m2 / z - m2Rec; 
    double r1Tmp    = sqrt(s1Tmp * s1Tmp + 4. * Q2 * m2Rec); 
    eNewRec         = 0.5 * (m2 + m2Rec + Q2) / m;
    pzNewRec        = -sideSign * 0.5 * r1Tmp / m;
    double pT2br    = Q2 * s3Tmp * (m2 / z - m2 - Q2) 
      - m2Sister * s1Tmp * s3Tmp - m2Rec * pow2(Q2 + m2Sister);
    if (pT2br <= 0.) return false; 
    pTbranch        = sqrt(pT2br) / r1Tmp;
    pzMother        = sideSign * (m * s3Tmp 
      - eNewRec * (m2 / z - Q2 - m2Rec - m2Sister)) / r1Tmp;
  }
  // Common final kinematics steps for both normal and rescattering.
  double eMother    = sqrt( pow2(pTbranch) + pow2(pzMother) );
  double pzSister   = pzMother + pzNewRec;
  double eSister    = sqrt( pow2(pTbranch) + pow2(pzSister) + m2Sister );
  Vec4 pMother( pTbranch, 0., pzMother, eMother );
  Vec4 pSister( pTbranch, 0., pzSister, eSister ); 
  Vec4 pNewRec(       0., 0., pzNewRec, eNewRec );

  // Current event and subsystem size.
  int eventSizeOld  = event.size();
  int systemSizeOld = partonSystemsPtr->sizeAll(iSysSel);

  // Save properties to be restored in case of user-hook veto of emission.
  int beamOff1 = 1 + beamOffset;
  int beamOff2 = 2 + beamOffset;
  int ev1Dau1V = event[beamOff1].daughter1();
  int ev2Dau1V = event[beamOff2].daughter1();
  vector<int> statusV, mother1V, mother2V, daughter1V, daughter2V;

  // Check if the first emission shoild be checked for removal
  bool canMergeFirst = (mergingHooksPtr != 0)
                     ? mergingHooksPtr->canVetoEmission() : false;
  if (canVetoEmission || canMergeFirst) {
    for ( int iCopy = 0; iCopy < systemSizeOld; ++iCopy) {
      int iOldCopy    = partonSystemsPtr->getAll(iSysSel, iCopy);
      statusV.push_back( event[iOldCopy].status());
      mother1V.push_back( event[iOldCopy].mother1());
      mother2V.push_back( event[iOldCopy].mother2());
      daughter1V.push_back( event[iOldCopy].daughter1());
      daughter2V.push_back( event[iOldCopy].daughter2());
    }  
  }

  // Take copy of existing system, to be given modified kinematics.
  // Incoming negative status. Rescattered also negative, but after copy.
  for ( int iCopy = 0; iCopy < systemSizeOld; ++iCopy) {
    int iOldCopy    = partonSystemsPtr->getAll(iSysSel, iCopy);
    int statusOld   = event[iOldCopy].status();
    int statusNew   = (iOldCopy == iDaughter 
      || iOldCopy == iRecoiler) ? statusOld : 44;
    int iNewCopy    = event.copy(iOldCopy, statusNew);
    if (statusOld < 0) event[iNewCopy].statusNeg();
  }
 
  // Define colour flow in branching. 
  // Default corresponds to f -> f + gamma.
  int colMother     = colDaughter;
  int acolMother    = acolDaughter;
  int colSister     = 0;
  int acolSister    = 0; 
  if (idSister == 22) ; 
  // q -> q + g and 50% of g -> g + g; need new colour.
  else if (idSister == 21 && ( (idMother > 0 && idMother < 9)
  || (idMother == 21 && rndmPtr->flat() < 0.5) ) ) {  
    colMother       = event.nextColTag();
    colSister       = colMother;
    acolSister      = colDaughter;
  // qbar -> qbar + g and other 50% of g -> g + g; need new colour.
  } else if (idSister == 21) {  
    acolMother      = event.nextColTag();
    acolSister      = acolMother;
    colSister       = acolDaughter;
  // q -> g + q.
  } else if (idDaughterNow == 21 && idMother > 0) { 
    colMother       = colDaughter;
    acolMother      = 0;
    colSister       = acolDaughter;
  // qbar -> g + qbar
  } else if (idDaughterNow == 21) {
    acolMother      = acolDaughter;
    colMother       = 0;
    acolSister      = colDaughter;
  // g -> q + qbar.
  } else if (idDaughterNow > 0 && idDaughterNow < 9) {
    acolMother      = event.nextColTag();
    acolSister      = acolMother;
  // g -> qbar + q.
  } else if (idDaughterNow < 0 && idDaughterNow > -9) {
    colMother       = event.nextColTag();
    colSister       = colMother;
  // q -> gamma + q.
  } else if (idDaughterNow == 22 && idMother > 0) {
    colMother       = event.nextColTag();
    colSister       = colMother; 
   // qbar -> gamma + qbar.
  } else if (idDaughterNow == 22) {
    acolMother      = event.nextColTag();
    acolSister      = acolMother;
  }   

  // Indices of partons involved. Add new sister.
  int iMother       = eventSizeOld + side - 1;
  int iNewRecoiler  = eventSizeOld + 2 - side;
  int iSister       = event.append( idSister, 43, iMother, 0, 0, 0,
     colSister, acolSister, pSister, mSister, sqrt(pT2) );

  // References to the partons involved.
  Particle& daughter    = event[iDaughter];
  Particle& mother      = event[iMother];
  Particle& newRecoiler = event[iNewRecoiler];
  Particle& sister      = event.back();

  // Replace old by new mother; update new recoiler.
  mother.id( idMother );
  mother.status( -41);
  mother.cols( colMother, acolMother);
  mother.p( pMother);
  newRecoiler.status( (normalRecoil) ? -42 : -46 );
  newRecoiler.p( pNewRec);
  if (!normalRecoil) newRecoiler.m( event[iRecoiler].m() ); 

  // Update mother and daughter pointers; also for beams.
  daughter.mothers( iMother, 0);
  mother.daughters( iSister, iDaughter); 
  if (iSysSel == 0) {
    event[beamOff1].daughter1( (side == 1) ? iMother : iNewRecoiler ); 
    event[beamOff2].daughter1( (side == 2) ? iMother : iNewRecoiler ); 
  }

  // Find boost to old rest frame.
  RotBstMatrix Mtot;
  if (normalRecoil) Mtot.bst(0., 0., (x2 - x1) / (x1 + x2) );
  else if (side == 1)
       Mtot.toCMframe( event[iDaughter].p(), event[iRecoiler].p() );
  else Mtot.toCMframe( event[iRecoiler].p(), event[iDaughter].p() );

  // Initially select phi angle of branching at random.
  double phi = 2. * M_PI * rndmPtr->flat();

  // Evaluate coefficient of azimuthal asymmetry from gluon polarization.
  findAsymPol( event, dipEndSel);
  int    iFinPol = dipEndSel->iFinPol;
  double cPol    = dipEndSel->asymPol;
  double phiPol  = (iFinPol == 0) ? 0. : event[iFinPol].phi(); 

  // If interference: try to match sister (anti)colour to final state.
  int    iFinInt = 0;
  double cInt    = 0.; 
  double phiInt  = 0.;
  if (doPhiIntAsym) {
    for (int i = 0; i < partonSystemsPtr->sizeOut(iSysSel); ++ i) {
      int iOut = partonSystemsPtr->getOut(iSysSel, i);
      if ( (acolSister != 0 && event[iOut].col() == acolSister) 
        || (colSister != 0 && event[iOut].acol() == colSister) ) 
        iFinInt = iOut;  
    }
    if (iFinInt != 0) {
      // Boost final-state parton to current frame of new kinematics.
      Vec4 pFinTmp = event[iFinInt].p();
      pFinTmp.rotbst(Mtot);
      double theFin = pFinTmp.theta();
      if (side == 2) theFin = M_PI - theFin;
      double theSis = pSister.theta();
      if (side == 2) theSis = M_PI - theSis;
      cInt = strengthIntAsym * 2. * theSis * theFin 
           / (pow2(theSis) + pow2(theFin));
      phiInt = event[iFinInt].phi();
    }
  }

  // Bias phi distribution for polarization and interference.
  if (iFinPol != 0 || iFinInt != 0) {
    double cPhiPol, cPhiInt, weight;
    do {
      phi     = 2. * M_PI * rndmPtr->flat();
      weight  = 1.;
      if (iFinPol !=0 ) {
        cPhiPol = cos(phi - phiPol);
        weight *= ( 1. + cPol * (2. * pow2(cPhiPol) - 1.) ) 
          / ( 1. + abs(cPol) );
      }
      if (iFinInt !=0 ) {
        cPhiInt = cos(phi - phiInt); 
        weight *= (1. - cInt) * (1. - cInt * cPhiInt)
          / (1. + pow2(cInt) - 2. * cInt * cPhiInt);
      }
    } while (weight < rndmPtr->flat());  
  }

  // Include rotation -phi on boost to old rest frame.
  Mtot.rot(0., -phi); 

  // Find boost from old rest frame to event cm frame.
  RotBstMatrix MfromRest;
  // The boost to the new rest frame.
  Vec4 sumNew       = pMother + pNewRec;
  double betaX      = sumNew.px() / sumNew.e();
  double betaZ      = sumNew.pz() / sumNew.e();
  MfromRest.bst( -betaX, 0., -betaZ);
  // Alignment of  radiator + recoiler to +- z axis, and rotation +phi.
  // Note: with spacelike (E < 0) recoiler p'_x_mother < 0 can happen!
  pMother.rotbst(MfromRest);  
  double theta = pMother.theta();
  if (pMother.px() < 0.) theta = -theta;
  if (side == 2) theta += M_PI;
  MfromRest.rot(-theta, phi); 
  // Boost to radiator + recoiler in event cm frame.
  if (normalRecoil) {
    MfromRest.bst( 0., 0., (x1New - x2New) / (x1New + x2New) );
  } else if (side == 1) {
    Vec4 pMotherWanted( 0., 0.,  0.5 * eCM * x1New, 0.5 * eCM * x1New);
    MfromRest.fromCMframe( pMotherWanted, event[iRecoiler].p() ); 

  } else {
    Vec4 pMotherWanted( 0., 0., -0.5 * eCM * x2New, 0.5 * eCM * x2New); 
    MfromRest.fromCMframe( event[iRecoiler].p(), pMotherWanted );
  }
  Mtot.rotbst(MfromRest);

  // Perform cumulative rotation/boost operation.
  // Mother, recoiler and sister from old rest frame to event cm frame.
  mother.rotbst(MfromRest);
  newRecoiler.rotbst(MfromRest);
  sister.rotbst(MfromRest);
  // The rest from (and to) event cm frame.
  for ( int i = eventSizeOld + 2; i < eventSizeOld + systemSizeOld; ++i) 
    event[i].rotbst(Mtot);  

  // Allow veto of branching. If so restore event record to before emission.
  if ( (canVetoEmission 
    && userHooksPtr->doVetoISREmission(eventSizeOld, event, iSysSel))
    || (canMergeFirst 
    && mergingHooksPtr->doVetoEmission( event )) ) {
    event.popBack( event.size() - eventSizeOld); 
    event[beamOff1].daughter1( ev1Dau1V);
    event[beamOff2].daughter1( ev2Dau1V);
    for ( int iCopy = 0; iCopy < systemSizeOld; ++iCopy) {
      int iOldCopy = partonSystemsPtr->getAll(iSysSel, iCopy);
      event[iOldCopy].status( statusV[iCopy]);
      event[iOldCopy].mothers( mother1V[iCopy], mother2V[iCopy]);
      event[iOldCopy].daughters( daughter1V[iCopy], daughter2V[iCopy]);
    }  
    return false;
  }

  // Update list of partons in system; adding newly produced one.
  partonSystemsPtr->setInA(iSysSel, eventSizeOld);
  partonSystemsPtr->setInB(iSysSel, eventSizeOld + 1);
  for (int iCopy = 2; iCopy < systemSizeOld; ++iCopy) 
    partonSystemsPtr->setOut(iSysSel, iCopy - 2, eventSizeOld + iCopy);
  partonSystemsPtr->addOut(iSysSel, eventSizeOld + systemSizeOld);
  partonSystemsPtr->setSHat(iSysSel, m2 / z);

  // Update info on radiating dipole ends (QCD or QED).
  for (int iDip = 0; iDip < int(dipEnd.size()); ++iDip)
  if ( dipEnd[iDip].system == iSysSel) {
    if (abs(dipEnd[iDip].side) == side) {    
      dipEnd[iDip].iRadiator = iMother;
      dipEnd[iDip].iRecoiler = iNewRecoiler;
      if (dipEnd[iDip].side > 0) 
        dipEnd[iDip].colType = mother.colType();
      else {
        dipEnd[iDip].chgType = 0;
        if ( (mother.isQuark() && doQEDshowerByQ)
          || (mother.isLepton() && doQEDshowerByL) ) 
          dipEnd[iDip].chgType = mother.chargeType();
      }
      // Kill ME corrections after first emission. 
      dipEnd[iDip].MEtype = 0;

    // Update info on recoiling dipole ends (QCD or QED).
    } else {
      dipEnd[iDip].iRadiator = iNewRecoiler;
      dipEnd[iDip].iRecoiler = iMother;
      // Optionally also kill recoiler ME corrections after first emission.
      if (!doMEafterFirst) dipEnd[iDip].MEtype = 0;
    }  
  }

  // Update info on beam remnants.
  BeamParticle& beamNow = (side == 1) ? *beamAPtr : *beamBPtr;
  double xNew = (side == 1) ? x1New : x2New;
  beamNow[iSysSel].update( iMother, idMother, xNew);
  // Redo choice of companion kind whenever new flavour.
  if (idMother != idDaughterNow) {
    beamNow.xfISR( iSysSel, idMother, xNew, factorMultFac * pT2);
    beamNow.pickValSeaComp();
  }
  BeamParticle& beamRec = (side == 1) ? *beamBPtr : *beamAPtr;
  beamRec[iSysSel].iPos( iNewRecoiler);

  // Store branching values of current dipole. (For rapidity ordering.)
  ++dipEndSel->nBranch;
  dipEndSel->pT2Old = pT2;
  dipEndSel->zOld   = z;

  // Update history if recoiler rescatters. 
  if (!normalRecoil) 
    event[iRecoilMother].daughters( iNewRecoiler, iNewRecoiler);

  // Start list of rescatterers that force changed kinematics.
  vector<int> iRescatterer;
  for ( int i = 0; i < systemSizeOld - 2; ++i) {
    int iOutNew = partonSystemsPtr->getOut( iSysSel, i);
    if (!event[iOutNew].isFinal()) iRescatterer.push_back(iOutNew);
  }

  // Start iterate over list of such rescatterers.
  int iRescNow = -1;
  while (++iRescNow < int(iRescatterer.size())) {

    // Identify partons that induce or are affected by rescatter shift.
    // In following Old is before change of kinematics, New after,
    // Out scatterer in outstate and In in instate of another system.
    // Daughter sequence is (iOutOld ->) iOutNew -> iInNew -> iInOld. 
    int iOutNew    = iRescatterer[iRescNow];    
    int iInOld     = event[iOutNew].daughter1();
    int iSysResc   = partonSystemsPtr->getSystemOf(iInOld, true);

    // Copy incoming partons of rescattered system and hook them up.
    int iOldA      = partonSystemsPtr->getInA(iSysResc);
    int iOldB      = partonSystemsPtr->getInB(iSysResc);
    bool rescSideA = event[iOldA].isRescatteredIncoming();
    int statusNewA = (rescSideA) ? -45 : -42;
    int statusNewB = (rescSideA) ? -42 : -45;
    int iNewA      = event.copy(iOldA, statusNewA);
    int iNewB      = event.copy(iOldB, statusNewB);

    // Copy outgoing partons of rescattered system and hook them up.
    int eventSize  = event.size();
    int sizeOutAB  = partonSystemsPtr->sizeOut(iSysResc);
    int iOldAB, statusOldAB, iNewAB;       
    for (int iOutAB = 0; iOutAB < sizeOutAB; ++iOutAB) {
      iOldAB       = partonSystemsPtr->getOut(iSysResc, iOutAB);
      statusOldAB  = event[iOldAB].status();
      iNewAB       = event.copy(iOldAB, 44);
      // Status could be negative for parton that rescatters in its turn.
      if (statusOldAB < 0) {
        event[iNewAB].statusNeg();     
        iRescatterer.push_back(iNewAB);
      }
    } 

    // Hook up new outgoing with new incoming parton.
    int iInNew     = (rescSideA) ? iNewA : iNewB;
    event[iOutNew].daughters( iInNew, iInNew);
    event[iInNew].mothers( iOutNew, iOutNew);

    // Rescale recoiling incoming parton for correct invariant mass.
    event[iInNew].p( event[iOutNew].p() );
    double momFac  = (rescSideA) 
                   ? event[iInOld].pPos() / event[iInNew].pPos()
                   : event[iInOld].pNeg() / event[iInNew].pNeg();
    int iInRec     = (rescSideA) ? iNewB : iNewA;

    // Rescatter: A previous boost may cause the light cone momentum of a
    //            rescattered parton to change sign. If this happens, tell
    //            parton level to try again.
    if (momFac < 0.0) {
      infoPtr->errorMsg("Warning in SpaceShower::branch: "
      "change in lightcone momentum sign; retrying parton level");
      rescatterFail = true;
      return false;
    }

    event[iInRec].rescale4( momFac);

    // Boost outgoing partons to new frame of incoming.
    RotBstMatrix MmodResc;
    MmodResc.toCMframe(  event[iOldA].p(), event[iOldB].p()); 
    MmodResc.fromCMframe(event[iNewA].p(), event[iNewB].p()); 
    for (int iOutAB = 0; iOutAB < sizeOutAB; ++iOutAB)
      event[eventSize + iOutAB].rotbst(MmodResc);

    // Update list of partons in system.
    partonSystemsPtr->setInA(iSysResc, iNewA);
    partonSystemsPtr->setInB(iSysResc, iNewB);
    for (int iCopy = 0; iCopy < sizeOutAB; ++iCopy) 
      partonSystemsPtr->setOut(iSysResc, iCopy, eventSize + iCopy);

    // Update info on radiating dipole ends (QCD or QED).
    for (int iDip = 0; iDip < int(dipEnd.size()); ++iDip)
    if ( dipEnd[iDip].system == iSysResc) {
      bool sideAnow = (abs(dipEnd[iDip].side) == 1); 
      dipEnd[iDip].iRadiator = (sideAnow) ? iNewA : iNewB;
      dipEnd[iDip].iRecoiler = (sideAnow) ? iNewB : iNewA;
    }

    // Update info on beam remnants.
    BeamParticle& beamResc = (rescSideA) ? *beamAPtr : *beamBPtr;
    beamResc[iSysResc].iPos( iInNew); 
    beamResc[iSysResc].p( event[iInNew].p() );
    beamResc[iSysResc].scaleX( 1. / momFac  );
    BeamParticle& beamReco = (rescSideA) ? *beamBPtr : *beamAPtr;
    beamReco[iSysResc].iPos( iInRec); 
    beamReco[iSysResc].scaleX( momFac);

  // End iterate over list of rescatterers.
  }

  // Check that beam momentum not used up by rescattered-system boosts.
    if (beamAPtr->xMax(-1) < 0.0 || beamBPtr->xMax(-1) < 0.0) {
      infoPtr->errorMsg("Warning in SpaceShower::branch: "
      "used up beam momentum; retrying parton level");
      rescatterFail = true;
      return false;
    }

  // Done without any errors.
  return true;

}

//--------------------------------------------------------------------------

// Find class of ME correction.

int SpaceShower::findMEtype( int iSys, Event& event) {

  // Default values and no action.
  int MEtype = 0; 
  if (!doMEcorrections) ;

  // Identify systems producing a single resonance.
  else if (partonSystemsPtr->sizeOut( iSys) == 1) {
    int idIn1 = event[partonSystemsPtr->getInA(iSys)].id();
    int idIn2 = event[partonSystemsPtr->getInA(iSys)].id();
    int idRes = event[partonSystemsPtr->getOut(iSys, 0)].id();
    if (iSys == 0) idResFirst  = abs(idRes);
    if (iSys == 1) idResSecond = abs(idRes);

    // f + fbar -> vector boson. 
    if ( (idRes == 23 || abs(idRes) == 24 || idRes == 32 
      || idRes == 33 || abs(idRes) == 34 || abs(idRes) == 41)
      && abs(idIn1) < 20 && abs(idIn2) < 20 ) MEtype = 1;

    // g + g, gamma + gamma  -> Higgs boson.
    if ( (idRes == 25 || idRes == 35 || idRes == 36) 
       && ( ( idIn1 == 21 && idIn2 == 21 ) 
       || ( idIn1 == 22 && idIn2 == 22 ) ) ) MEtype = 2; 

    // f + fbar  -> Higgs boson.
    if ( (idRes == 25 || idRes == 35 || idRes == 36) 
       && abs(idIn1) < 20 && abs(idIn2) < 20 ) MEtype = 3; 
  }

  // Done.
  return MEtype;

}

//--------------------------------------------------------------------------

// Provide maximum of expected ME weight; for preweighting of evolution.

double SpaceShower::calcMEmax( int MEtype, int idMother, int idDaughterIn) {

  // Currently only one non-unity case: g(gamma) f -> V f'.
  if (MEtype == 1 && idMother > 20 && idDaughterIn < 20) return 3.;
  return 1.;

}  

//--------------------------------------------------------------------------

// Provide actual ME weight for current branching.
// Note: currently ME corrections are only allowed for first branching 
// on each side, so idDaughter is essentially known and checks overkill.

double SpaceShower::calcMEcorr(int MEtype, int idMother, int idDaughterIn,
  double M2, double z, double Q2) {

  // Convert to Mandelstam variables. Sometimes may need to swap later.
  double sH = M2 / z;
  double tH = -Q2;
  double uH = Q2 - M2 * (1. - z) / z;
  int idMabs = abs(idMother);
  int idDabs = abs(idDaughterIn);

  // Corrections for f + fbar -> s-channel vector boson.
  if (MEtype == 1) {
    if (idMabs < 20 && idDabs < 20) {
      return (tH*tH + uH*uH + 2. * M2 * sH) / (sH*sH + M2*M2); 
    } else if (idDabs < 20) {
      // g(gamma) f -> V f': -Q2 = (p_g - p_f')^2 in PS while 
      // tHat = (p_f - p_f')^2 in ME so need to swap tHat <-> uHat.  
      swap( tH, uH); 
      return (sH*sH + uH*uH + 2. * M2 * tH) / (pow2(sH - M2) + M2*M2); 
    }

  // Corrections for g + g -> Higgs boson.
  } else if (MEtype == 2) {
    if (idMabs < 20 && idDabs > 20) {
      return (sH*sH + uH*uH) / (sH*sH + pow2(sH - M2)); 
    } else if (idDabs > 20) {
      return 0.5 * (pow4(sH) + pow4(tH) + pow4(uH) + pow4(M2)) 
        / pow2(sH*sH - M2 * (sH - M2)); 
    }    

  // Corrections for f + fbar -> Higgs boson (f = b mainly).
  } else if (MEtype == 3) {
    if (idMabs < 20 && idDabs < 20) {
      // The PS and ME answers agree for f fbar -> H g/gamma. 
      return 1.; 
    } else if (idDabs < 20) {
      // Need to swap tHat <-> uHat, cf. vector-boson production above. 
      swap( tH, uH); 
      return (sH*sH + uH*uH + 2. * (M2 - uH) * (M2 - sH)) 
             / (pow2(sH - M2) + M2*M2); 
    }    
  }

  return 1.;

}

//--------------------------------------------------------------------------

// Find coefficient of azimuthal asymmetry from gluon polarization.

void SpaceShower::findAsymPol( Event& event, SpaceDipoleEnd* dip) {

  // Default is no asymmetry. Only gluons are studied.
  dip->iFinPol   = 0;
  dip->asymPol   = 0.;
  int iRad       = dip->iRadiator;
  if (!doPhiPolAsym || dip->idDaughter != 21) return;

  // At least two particles in final state, whereof at least one coloured.
  int systemSizeOut = partonSystemsPtr->sizeOut( iSysSel);
  if (systemSizeOut < 2) return;
  bool foundColOut  = false;
  for (int ii = 0; ii < systemSizeOut; ++ii) {
    int i = partonSystemsPtr->getOut( iSysSel, ii); 
    if (event[i].col() != 0 || event[i].acol() != 0) foundColOut = true;
  }
  if (!foundColOut) return;

  // Check if granddaughter in final state of hard scattering.
  // (May need to trace across carbon copies to find granddaughters.)
  // If so, only accept 2 -> 2 scatterings with gg or qq in final state.
  int iGrandD1 = event[iRad].daughter1();
  int iGrandD2 = event[iRad].daughter2();
  bool traceCopy = false;
  do {
    traceCopy = false;
    if (iGrandD1 > 0 && iGrandD2 == iGrandD1) {
      iGrandD1 = event[iGrandD2].daughter1();
      iGrandD2 = event[iGrandD2].daughter2();
      traceCopy = true;
    }
  } while (traceCopy);
  int statusGrandD1 = event[ iGrandD1 ].statusAbs();
  bool isHardProc  = (statusGrandD1 == 23 || statusGrandD1 == 33); 
  if (isHardProc) {
    if (iGrandD2 != iGrandD1 + 1) return; 
    if (event[iGrandD1].isGluon() && event[iGrandD2].isGluon());
    else if (event[iGrandD1].isQuark() && event[iGrandD2].isQuark());
    else return;
  }
  dip->iFinPol = iGrandD1;   

  // Coefficient from gluon production.
  if (dip->idMother == 21) dip->asymPol = pow2( (1. - dip->z) 
    / (1. - dip->z * (1. - dip->z) ) );
  else dip->asymPol = 2. * (1. - dip->z) / (1. + pow2(1. - dip->z) );

  // Coefficients from gluon decay. Put z = 1/2 for hard process.
  double zDau  = (isHardProc) ? 0.5 : dip->zOld;
  if (event[iGrandD1].isGluon()) dip->asymPol *= pow2( (1. - zDau) 
    / (1. - zDau * (1. - zDau) ) );
  else  dip->asymPol *= -2. * zDau *( 1. - zDau ) 
    / (1. - 2. * zDau * (1. - zDau) );

}

//--------------------------------------------------------------------------

// Print the list of dipoles.

void SpaceShower::list(ostream& os) const {

  // Header.
  os << "\n --------  PYTHIA SpaceShower Dipole Listing  -------------- \n"
     << "\n    i  syst  side   rad   rec       pTmax  col  chg   ME rec \n"
     << fixed << setprecision(3);
  
  // Loop over dipole list and print it.
  for (int i = 0; i < int(dipEnd.size()); ++i) 
  os << setw(5) << i << setw(6) << dipEnd[i].system 
     << setw(6) << dipEnd[i].side << setw(6) << dipEnd[i].iRadiator 
     << setw(6) << dipEnd[i].iRecoiler << setw(12) << dipEnd[i].pTmax 
     << setw(5) << dipEnd[i].colType << setw(5) << dipEnd[i].chgType
     << setw(5) << dipEnd[i].MEtype << setw(4) 
     << dipEnd[i].normalRecoil << "\n";

  // Done.
  os << "\n --------  End PYTHIA SpaceShower Dipole Listing  ----------"
     << endl;
  

}
 
//==========================================================================

} // end namespace Pythia8

