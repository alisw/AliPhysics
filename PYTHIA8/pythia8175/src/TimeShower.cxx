// TimeShower.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the TimeShower class. 

#include "TimeShower.h"

namespace Pythia8 {

//==========================================================================

// The TimeShower class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// For small x approximate 1 - sqrt(1 - x) by x/2.
const double TimeShower::SIMPLIFYROOT = 1e-8;

// Do not allow x too close to 0 or 1 in matrix element expressions.
// Warning: cuts into phase space for E_CM > 2 * pTmin * sqrt(1/XMARGIN),
// i.e. will become problem roughly for E_CM > 10^6 GeV.
const double TimeShower::XMARGIN      = 1e-12;
const double TimeShower::XMARGINCOMB  = 1e-4;

// Lower limit on PDF value in order to avoid division by zero.
const double TimeShower::TINYPDF      = 1e-10;

// Big starting value in search for smallest invariant-mass pair.
const double TimeShower::LARGEM2      = 1e20;

// In g -> q qbar or gamma -> f fbar require m2_pair > this * m2_q/f. 
const double TimeShower::THRESHM2      = 4.004;

// Never pick pT so low that alphaS is evaluated too close to Lambda_3. 
const double TimeShower::LAMBDA3MARGIN = 1.1;

// Rescatter: rescattering + ISR + FSR + primordial kT can lead to
//            systems not locally conserving momentum.
// Fix up momentum in intermediate systems with rescattering
const bool   TimeShower::FIXRESCATTER          = true;
// Veto negative energies when using FIXRESCATTER option.
const bool   TimeShower::VETONEGENERGY         = false;
// Do not allow too large time- or spacelike virtualities in fixing-up.
const double TimeShower::MAXVIRTUALITYFRACTION = 0.5;
// Do not allow too large negative spacelike energy in system rest frame.
const double TimeShower::MAXNEGENERGYFRACTION  = 0.7; 

//--------------------------------------------------------------------------

// Initialize alphaStrong, alphaEM and related pTmin parameters.

void TimeShower::init( BeamParticle* beamAPtrIn, 
  BeamParticle* beamBPtrIn) {

  // Store input pointers for future use. 
  beamAPtr           = beamAPtrIn;
  beamBPtr           = beamBPtrIn;

  // Main flags.
  doQCDshower        = settingsPtr->flag("TimeShower:QCDshower");
  doQEDshowerByQ     = settingsPtr->flag("TimeShower:QEDshowerByQ");
  doQEDshowerByL     = settingsPtr->flag("TimeShower:QEDshowerByL");
  doQEDshowerByGamma = settingsPtr->flag("TimeShower:QEDshowerByGamma");
  doMEcorrections    = settingsPtr->flag("TimeShower:MEcorrections");
  doMEafterFirst     = settingsPtr->flag("TimeShower:MEafterFirst");
  doPhiPolAsym       = settingsPtr->flag("TimeShower:phiPolAsym"); 
  doInterleave       = settingsPtr->flag("TimeShower:interleave"); 
  allowBeamRecoil    = settingsPtr->flag("TimeShower:allowBeamRecoil"); 
  dampenBeamRecoil   = settingsPtr->flag("TimeShower:dampenBeamRecoil"); 
  recoilToColoured   = settingsPtr->flag("TimeShower:recoilToColoured"); 

  // Matching in pT of hard interaction or MPI to shower evolution.
  pTmaxMatch         = settingsPtr->mode("TimeShower:pTmaxMatch"); 
  pTdampMatch        = settingsPtr->mode("TimeShower:pTdampMatch"); 
  pTmaxFudge         = settingsPtr->parm("TimeShower:pTmaxFudge"); 
  pTmaxFudgeMPI      = settingsPtr->parm("TimeShower:pTmaxFudgeMPI"); 
  pTdampFudge        = settingsPtr->parm("TimeShower:pTdampFudge"); 

  // Charm and bottom mass thresholds.
  mc                 = particleDataPtr->m0(4); 
  mb                 = particleDataPtr->m0(5); 
  m2c                = mc * mc;
  m2b                = mb * mb;

  // Parameters of scale choices.
  renormMultFac     = settingsPtr->parm("TimeShower:renormMultFac");
  factorMultFac     = settingsPtr->parm("TimeShower:factorMultFac");
       
  // Parameters of alphaStrong generation.
  alphaSvalue        = settingsPtr->parm("TimeShower:alphaSvalue");
  alphaSorder        = settingsPtr->mode("TimeShower:alphaSorder");
  alphaS2pi          = 0.5 * alphaSvalue / M_PI;

  // Initialize alphaStrong generation.
  alphaS.init( alphaSvalue, alphaSorder); 
  
  // Lambda for 5, 4 and 3 flavours.
  Lambda3flav        = alphaS.Lambda3(); 
  Lambda4flav        = alphaS.Lambda4(); 
  Lambda5flav        = alphaS.Lambda5(); 
  Lambda5flav2       = pow2(Lambda5flav);
  Lambda4flav2       = pow2(Lambda4flav);
  Lambda3flav2       = pow2(Lambda3flav);
 
  // Parameters of QCD evolution. Warn if pTmin must be raised.
  nGluonToQuark      = settingsPtr->mode("TimeShower:nGluonToQuark");
  pTcolCutMin        = settingsPtr->parm("TimeShower:pTmin");
  if (pTcolCutMin > LAMBDA3MARGIN * Lambda3flav / sqrt(renormMultFac)) 
    pTcolCut         = pTcolCutMin;
  else { 
    pTcolCut         = LAMBDA3MARGIN * Lambda3flav / sqrt(renormMultFac);
    ostringstream newPTcolCut;
    newPTcolCut << fixed << setprecision(3) << pTcolCut;
    infoPtr->errorMsg("Warning in TimeShower::init: pTmin too low",
	 	      ", raised to " + newPTcolCut.str() );
    infoPtr->setTooLowPTmin(true);
  }
  pT2colCut          = pow2(pTcolCut);  
       
  // Parameters of alphaEM generation .
  alphaEMorder       = settingsPtr->mode("TimeShower:alphaEMorder");

  // Initialize alphaEM generation.
  alphaEM.init( alphaEMorder, settingsPtr); 
 
  // Parameters of QED evolution.
  nGammaToQuark      = settingsPtr->mode("TimeShower:nGammaToQuark");
  nGammaToLepton     = settingsPtr->mode("TimeShower:nGammaToLepton");
  pTchgQCut          = settingsPtr->parm("TimeShower:pTminChgQ"); 
  pT2chgQCut         = pow2(pTchgQCut);
  pTchgLCut          = settingsPtr->parm("TimeShower:pTminChgL"); 
  pT2chgLCut         = pow2(pTchgLCut);
  mMaxGamma          = settingsPtr->parm("TimeShower:mMaxGamma"); 
  m2MaxGamma         = pow2(mMaxGamma);

  // Consisteny check for gamma -> f fbar variables.
  if (nGammaToQuark <= 0 && nGammaToLepton <= 0) doQEDshowerByGamma = false; 

  // Possibility of a global recoil stategy, e.g. for MC@NLO.
  globalRecoil       = settingsPtr->flag("TimeShower:globalRecoil");
  nMaxGlobalRecoil   = settingsPtr->mode("TimeShower:nMaxGlobalRecoil");

  // Fraction and colour factor of gluon emission off onium octat state.
  octetOniumFraction = settingsPtr->parm("TimeShower:octetOniumFraction");
  octetOniumColFac   = settingsPtr->parm("TimeShower:octetOniumColFac");

  // Z0 properties needed for gamma/Z0 mixing.
  mZ                 = particleDataPtr->m0(23);
  gammaZ             = particleDataPtr->mWidth(23);
  thetaWRat          = 1. / (16. * coupSMPtr->sin2thetaW() 
                       * coupSMPtr->cos2thetaW());

  // May have to fix up recoils related to rescattering. 
  allowRescatter     = settingsPtr->flag("PartonLevel:MPI") 
    && settingsPtr->flag("MultipartonInteractions:allowRescatter");

  // Hidden Valley scenario with further shower activity.
  doHVshower         = settingsPtr->flag("HiddenValley:FSR");
  nCHV               = settingsPtr->mode("HiddenValley:Ngauge");
  alphaHVfix         = settingsPtr->parm("HiddenValley:alphaFSR");
  pThvCut            = settingsPtr->parm("HiddenValley:pTminFSR");
  pT2hvCut           = pThvCut * pThvCut; 
  CFHV               = (nCHV == 1) ? 1. : (nCHV * nCHV - 1.)/(2. * nCHV); 
  idHV               = (nCHV == 1) ? 4900022 : 4900021;
  mHV                = particleDataPtr->m0(idHV);
  brokenHVsym        = (nCHV == 1 && mHV > 0.);

  // Possibility of two predetermined hard emissions in event.
  doSecondHard    = settingsPtr->flag("SecondHard:generate");

  // Possibility to allow user veto of emission step.
  canVetoEmission    = (userHooksPtr != 0) 
                     ? userHooksPtr->canVetoFSREmission() : false;

}

//--------------------------------------------------------------------------

// Find whether to limit maximum scale of emissions.
// Also allow for dampening at factorization or renormalization scale. 

bool TimeShower::limitPTmax( Event& event, double Q2Fac, double Q2Ren) {

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

// Top-level routine to do a full time-like shower in resonance decay.

int TimeShower::shower( int iBeg, int iEnd, Event& event, double pTmax,
  int nBranchMax) {

  // Add new system, automatically with two empty beam slots.
  int iSys = partonSystemsPtr->addSys();
    
  // Loop over allowed range to find all final-state particles.
  Vec4 pSum;
  for (int i = iBeg; i <= iEnd; ++i) if (event[i].isFinal()) {
    partonSystemsPtr->addOut( iSys, i);
    pSum += event[i].p();
  }
  partonSystemsPtr->setSHat( iSys, pSum.m2Calc() );

  // Let prepare routine do the setup.    
  prepare( iSys, event, true);

  // Begin evolution down in pT from hard pT scale. 
  int nBranch  = 0;
  pTLastBranch = 0.;
  do {
    double pTtimes = pTnext( event, pTmax, 0.);

    // Do a final-state emission (if allowed).
    if (pTtimes > 0.) {
      if (branch( event)) {
        ++nBranch;
        pTLastBranch = pTtimes;
      }
      pTmax = pTtimes;
    }
    
    // Keep on evolving until nothing is left to be done.
    else pTmax = 0.;
  } while (pTmax > 0. && (nBranchMax <= 0 || nBranch < nBranchMax)); 

  // Return number of emissions that were performed.
  return nBranch;

}

//--------------------------------------------------------------------------

// Top-level routine for QED radiation in hadronic decay to two leptons.
// Intentionally only does photon radiation, i.e. no photon branchings. 

int TimeShower::showerQED( int i1, int i2, Event& event, double pTmax) {

  // Add new system, automatically with two empty beam slots.
  int iSys = partonSystemsPtr->addSys();
  partonSystemsPtr->addOut( iSys, i1);
  partonSystemsPtr->addOut( iSys, i2);
  partonSystemsPtr->setSHat( iSys, m2(event[i1], event[i2]) );

  // Charge type of two leptons tells whether MEtype is gamma*/Z0 or W+-.
  int iChg1  = event[i1].chargeType();
  int iChg2  = event[i2].chargeType();
  int MEtype = (iChg1 + iChg2 == 0) ? 102 : 101;

  // Fill dipole-ends list.    
  dipEnd.resize(0);
  if (iChg1 != 0) dipEnd.push_back( TimeDipoleEnd(i1, i2, pTmax,
      0, iChg1, 0, 0, iSys, MEtype, i2) );
  if (iChg2 != 0) dipEnd.push_back( TimeDipoleEnd(i2, i1, pTmax,
      0, iChg2, 0, 0, iSys, MEtype, i1) );

  // Begin evolution down in pT from hard pT scale. 
  int nBranch  = 0;
  pTLastBranch = 0.;
  do {

    // Begin loop over all possible radiating dipole ends.
    dipSel  = 0;
    iDipSel = -1;
    double pT2sel = 0.;
    for (int iDip = 0; iDip < int(dipEnd.size()); ++iDip) {
      TimeDipoleEnd& dip = dipEnd[iDip]; 

      // Dipole properties. 
      dip.mRad  = event[dip.iRadiator].m(); 
      dip.mRec  = event[dip.iRecoiler].m(); 
      dip.mDip  = m( event[dip.iRadiator], event[dip.iRecoiler] );
      dip.m2Rad = pow2(dip.mRad);
      dip.m2Rec = pow2(dip.mRec);
      dip.m2Dip = pow2(dip.mDip);

      // Find maximum evolution scale for dipole.
      dip.m2DipCorr    = pow2(dip.mDip - dip.mRec) - dip.m2Rad; 
      double pTbegDip  = min( pTmax, dip.pTmax ); 
      double pT2begDip = min( pow2(pTbegDip), 0.25 * dip.m2DipCorr);

      // Do QED evolution where relevant.
      dip.pT2 = 0.;
      if (pT2begDip > pT2sel) {
        pT2nextQED( pT2begDip, pT2sel, dip, event);

        // Update if found larger pT than current maximum. End dipole loop.
        if (dip.pT2 > pT2sel) {
          pT2sel  = dip.pT2;
          dipSel  = &dip;
          iDipSel = iDip;
        }
      } 
    }
    double pTsel = (dipSel == 0) ? 0. : sqrt(pT2sel); 

    // Do a final-state emission (if allowed).
    if (pTsel > 0.) {

      // Find initial radiator and recoiler particles in dipole branching.
      int iRadBef      = dipSel->iRadiator;
      int iRecBef      = dipSel->iRecoiler;
      Particle& radBef = event[iRadBef]; 
      Particle& recBef = event[iRecBef];
      Vec4 pRadBef     = event[iRadBef].p(); 
      Vec4 pRecBef     = event[iRecBef].p();  
      
      // Construct kinematics in dipole rest frame; massless emitter.
      double pTorig       = sqrt( dipSel->pT2);
      double eRadPlusEmt  = 0.5 * (dipSel->m2Dip + dipSel->m2 - dipSel->m2Rec) 
        / dipSel->mDip;
      double e2RadPlusEmt = pow2(eRadPlusEmt);
      double pzRadPlusEmt = 0.5 * sqrtpos( pow2(dipSel->m2Dip - dipSel->m2 
        - dipSel->m2Rec) - 4. * dipSel->m2 * dipSel->m2Rec ) / dipSel->mDip;
      double pT2corr = dipSel->m2 * (e2RadPlusEmt * dipSel->z 
        * (1. - dipSel->z) - 0.25 * dipSel->m2) / pow2(pzRadPlusEmt);
      double pTcorr       = sqrtpos( pT2corr );
      double pzRad        = (e2RadPlusEmt * dipSel->z - 0.5 * dipSel->m2) 
        / pzRadPlusEmt;
      double pzEmt        = (e2RadPlusEmt * (1. - dipSel->z) 
        - 0.5 * dipSel->m2) / pzRadPlusEmt;
      double mRad         = dipSel->mRad;
      double mEmt         = 0.;

      // Kinematics reduction for radiator mass. 
      double m2Ratio    = dipSel->m2Rad / dipSel->m2;
      pTorig           *= 1. - m2Ratio; 
      pTcorr           *= 1. - m2Ratio; 
      pzRad            += pzEmt * m2Ratio;
      pzEmt            *= 1. - m2Ratio; 

      // Store kinematics of branching in dipole rest frame.
      double phi = 2. * M_PI * rndmPtr->flat();
      Vec4 pRad = Vec4( pTcorr * cos(phi), pTcorr * sin(phi), pzRad, 
        sqrt( pow2(pTcorr) + pow2(pzRad) + pow2(mRad) ) );
      Vec4 pEmt = Vec4( -pRad.px(), -pRad.py(), pzEmt,
        sqrt( pow2(pTcorr) + pow2(pzEmt) + pow2(mEmt) ) );
      Vec4 pRec = Vec4( 0., 0., -pzRadPlusEmt, 
        sqrt( pow2(pzRadPlusEmt) + dipSel->m2Rec ) );

      // Rotate and boost dipole products to the event frame.
      RotBstMatrix M;
      M.fromCMframe(pRadBef, pRecBef);
      pRad.rotbst(M);
      pEmt.rotbst(M);
      pRec.rotbst(M);

      // Define new particles from dipole branching.
      Particle rad = Particle(radBef.id(), 51, iRadBef, 0, 0, 0, 
        radBef.col(), radBef.acol(), pRad, mRad, pTsel); 
      Particle emt = Particle(22, 51, iRadBef, 0, 0, 0,
        0, 0, pEmt, mEmt, pTsel);
      Particle rec = Particle(recBef.id(),  52, iRecBef, iRecBef, 0, 0, 
        recBef.col(), recBef.acol(), pRec, dipSel->mRec, pTsel); 

      // ME corrections can lead to branching being rejected.
      if (dipSel->MEtype == 0 
        || findMEcorr( dipSel, rad, rec, emt, false) > rndmPtr->flat() ) { 
 
        // Shower may occur at a displaced vertex, or for unstable particle.
        if (radBef.hasVertex()) {
          rad.vProd( radBef.vProd() );
          emt.vProd( radBef.vProd() );
        }
        if (recBef.hasVertex()) rec.vProd( recBef.vProd() );
        rad.tau( event[iRadBef].tau() );
        rec.tau( event[iRecBef].tau() );

        // Put new particles into the event record.
        int iRad = event.append(rad);
        int iEmt = event.append(emt);
        event[iRadBef].statusNeg();
        event[iRadBef].daughters( iRad, iEmt); 
        int iRec = event.append(rec);
        event[iRecBef].statusNeg();
        event[iRecBef].daughters( iRec, iRec);

        // Update to new dipole ends.
        dipSel->iRadiator = iRad;
        dipSel->iRecoiler = iRec;
        dipSel->pTmax = pTsel;

        // Update other dipoles that also involved the radiator or recoiler.
        for (int i = 0; i < int(dipEnd.size()); ++i) if (i != iDipSel) {
          if (dipEnd[i].iRadiator  == iRadBef) dipEnd[i].iRadiator  = iRad;
          if (dipEnd[i].iRecoiler  == iRadBef) dipEnd[i].iRecoiler  = iRad;
          if (dipEnd[i].iMEpartner == iRadBef) dipEnd[i].iMEpartner = iRad;
          if (dipEnd[i].iRadiator  == iRecBef) dipEnd[i].iRadiator  = iRec;
          if (dipEnd[i].iRecoiler  == iRecBef) dipEnd[i].iRecoiler  = iRec;
          if (dipEnd[i].iMEpartner == iRecBef) dipEnd[i].iMEpartner = iRec;
        }

        // Done with branching      
        ++nBranch;
        pTLastBranch = pTsel;
      }
      pTmax = pTsel;
    }
    
    // Keep on evolving until nothing is left to be done.
    else pTmax = 0.;
  } while (pTmax > 0.); 

  // Return number of emissions that were performed.
  return nBranch;

}

//--------------------------------------------------------------------------

// Prepare system for evolution; identify ME.

void TimeShower::prepare( int iSys, Event& event, bool limitPTmaxIn) {

  // Reset dipole-ends list for first interaction and for resonance decays.
  int iInA = partonSystemsPtr->getInA(iSys);
  int iInB = partonSystemsPtr->getInB(iSys); 
  if (iSys == 0 || iInA == 0) dipEnd.resize(0);
  int dipEndSizeBeg = dipEnd.size();

  // No dipoles for 2 -> 1 processes.
  if (partonSystemsPtr->sizeOut(iSys) < 2) return;
 
  // In case of DPS overwrite limitPTmaxIn by saved value.
  if (doSecondHard && iSys == 0) limitPTmaxIn = dopTlimit1; 
  if (doSecondHard && iSys == 1) limitPTmaxIn = dopTlimit2; 

  // Loop through final state of system to find possible dipole ends.
  for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
    int iRad = partonSystemsPtr->getOut( iSys, i);
    if (event[iRad].isFinal() && event[iRad].scale() > 0.) {

      // Identify colour octet onium state. Check whether QCD shower allowed.
      int idRad    = event[iRad].id();
      int idRadAbs = abs(idRad);
      bool isOctetOnium 
        = ( idRad == 9900441 || idRad == 9900443 || idRad == 9910441 
         || idRad == 9900551 || idRad == 9900553 || idRad == 9910551 ); 
      bool doQCD = doQCDshower;
      if (doQCD && isOctetOnium) 
        doQCD = (rndmPtr->flat() < octetOniumFraction);

      // Find dipole end formed by colour index.
      int colTag = event[iRad].col();    
      if (doQCD && colTag > 0) setupQCDdip( iSys, i,  colTag,  1, event, 
        isOctetOnium, limitPTmaxIn); 

      // Find dipole end formed by anticolour index.
      int acolTag = event[iRad].acol();     
      if (doQCD && acolTag > 0) setupQCDdip( iSys, i, acolTag, -1, event, 
        isOctetOnium, limitPTmaxIn); 

      // Find "charge-dipole" and "photon-dipole" ends. 
      int  chgType  = event[iRad].chargeType();  
      bool doChgDip = (chgType != 0) 
                       && ( ( doQEDshowerByQ && event[iRad].isQuark()  )
                         || ( doQEDshowerByL && event[iRad].isLepton() ) );
      int  gamType  = (idRad == 22) ? 1 : 0;
      bool doGamDip = (gamType == 1) && doQEDshowerByGamma;
      if (doChgDip || doGamDip) setupQEDdip( iSys, i, chgType, gamType, 
	 event, limitPTmaxIn);

      // Find Hidden Valley dipole ends.
      bool isHVrad =  (idRadAbs > 4900000 && idRadAbs < 4900007)
                   || (idRadAbs > 4900010 && idRadAbs < 4900017)
                   || idRadAbs == 4900101;  
      if (doHVshower && isHVrad) setupHVdip( iSys, i, event, limitPTmaxIn); 

    // End loop over system final state. Have now found the dipole ends.
    }
  }

  // Loop through dipole ends to find matrix element corrections.
  for (int iDip = dipEndSizeBeg; iDip < int(dipEnd.size()); ++iDip) 
    findMEtype( event, dipEnd[iDip]); 

  // Update dipole list after a multiparton interactions rescattering.
  if (iSys > 0 && ( (iInA > 0 && event[iInA].status() == -34) 
    || (iInB > 0 && event[iInB].status() == -34) ) ) 
    rescatterUpdate( iSys, event); 

}

//--------------------------------------------------------------------------

// Update dipole list after a multiparton interactions rescattering. 
void TimeShower::rescatterUpdate( int iSys, Event& event) {
 
  // Loop over two incoming partons in system; find their rescattering mother.
  // (iOut is outgoing from old system = incoming iIn of rescattering system.) 
  for (int iResc = 0; iResc < 2; ++iResc) { 
    int iIn = (iResc == 0) ? partonSystemsPtr->getInA(iSys)
                           : partonSystemsPtr->getInB(iSys); 
    if (iIn == 0 || event[iIn].status() != -34) continue; 
    int iOut = event[iIn].mother1();

    // Loop over all dipoles.
    int dipEndSize = dipEnd.size();
    for (int iDip = 0; iDip < dipEndSize; ++iDip) {
      TimeDipoleEnd& dipNow = dipEnd[iDip];

      // Kill dipoles where rescattered parton is radiator.
      if (dipNow.iRadiator == iOut) {
        dipNow.colType = 0;
        dipNow.chgType = 0;
        dipNow.gamType = 0;
        continue;
      }
      // No matrix element for dipoles between scatterings.
      if (dipNow.iMEpartner == iOut) {
        dipNow.MEtype     =  0;
        dipNow.iMEpartner = -1;
      }

      // Update dipoles where outgoing rescattered parton is recoiler.
      if (dipNow.iRecoiler == iOut) {
        int iRad = dipNow.iRadiator;

        // Colour dipole: recoil in final state, initial state or new.
        if (dipNow.colType > 0) {
          int colTag = event[iRad].col(); 
          bool done  = false;
          for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
            int iRecNow = partonSystemsPtr->getOut( iSys, i);  
            if (event[iRecNow].acol() == colTag) {
              dipNow.iRecoiler = iRecNow;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              done             = true;
              break;
            }
          }
          if (!done) {
            int iIn2 = (iResc == 0) ? partonSystemsPtr->getInB(iSys)
                                    : partonSystemsPtr->getInA(iSys); 
            if (event[iIn2].col() == colTag) {
              dipNow.iRecoiler = iIn2;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              int isrType      = event[iIn2].mother1();
              // This line in case mother is a rescattered parton.
              while (isrType > 2 + beamOffset) 
                isrType = event[isrType].mother1();
              if (isrType > 2) isrType -= beamOffset;
              dipNow.isrType   = isrType;
              done             = true;
            }
          }
          // If above options failed, then create new dipole.
          if (!done) {
            int iRadNow = partonSystemsPtr->getIndexOfOut(dipNow.system, iRad);
            if (iRadNow != -1)
              setupQCDdip(dipNow.system, iRadNow, event[iRad].col(), 1,
                          event, dipNow.isOctetOnium, true);
            else
              infoPtr->errorMsg("Warning in TimeShower::rescatterUpdate: "
              "failed to locate radiator in system");

            dipNow.colType = 0;
            dipNow.chgType = 0;
            dipNow.gamType = 0; 

            infoPtr->errorMsg("Warning in TimeShower::rescatterUpdate: "
            "failed to locate new recoiling colour partner");
          }

        // Anticolour dipole: recoil in final state, initial state or new.
        } else if (dipNow.colType < 0) {
          int  acolTag = event[iRad].acol(); 
          bool done    = false;  
          for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
            int iRecNow = partonSystemsPtr->getOut( iSys, i);  
            if (event[iRecNow].col() == acolTag) {
              dipNow.iRecoiler = iRecNow;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              done             = true;
              break;
            }
          }
          if (!done) {
            int iIn2 = (iResc == 0) ? partonSystemsPtr->getInB(iSys)
                                    : partonSystemsPtr->getInA(iSys); 
            if (event[iIn2].acol() == acolTag) {
              dipNow.iRecoiler = iIn2;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              int isrType      = event[iIn2].mother1();
              // This line in case mother is a rescattered parton.
              while (isrType > 2 + beamOffset) 
                isrType = event[isrType].mother1();
              if (isrType > 2) isrType -= beamOffset;
              dipNow.isrType   = isrType;
              done             = true;
            }
          } 
          // If above options failed, then create new dipole.
          if (!done) {
            int iRadNow = partonSystemsPtr->getIndexOfOut(dipNow.system, iRad);
            if (iRadNow != -1)
              setupQCDdip(dipNow.system, iRadNow, event[iRad].acol(), -1,
                          event, dipNow.isOctetOnium, true);
            else
              infoPtr->errorMsg("Warning in TimeShower::rescatterUpdate: "
              "failed to locate radiator in system");

            dipNow.colType = 0;
            dipNow.chgType = 0;
            dipNow.gamType = 0;

            infoPtr->errorMsg("Warning in TimeShower::rescatterUpdate: "
            "failed to locate new recoiling colour partner");
          }

        // Charge or photon dipoles: same flavour in final or initial state.
        } else if (dipNow.chgType != 0 || dipNow.gamType != 0) {
          int  idTag = event[dipNow.iRecoiler].id();
          bool done  = false;  
          for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
            int iRecNow = partonSystemsPtr->getOut( iSys, i); 
            if (event[iRecNow].id() == idTag) {
              dipNow.iRecoiler = iRecNow;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              done             = true;
              break;
            }
          }
          if (!done) {
            int iIn2 = (iResc == 0) ? partonSystemsPtr->getInB(iSys)
                                    : partonSystemsPtr->getInA(iSys); 
            if (event[iIn2].id() == -idTag) {
              dipNow.iRecoiler = iIn2;
              dipNow.systemRec = iSys;
              dipNow.MEtype    = 0;
              int isrType      = event[iIn2].mother1();
              // This line in case mother is a rescattered parton.
              while (isrType > 2 + beamOffset) 
                isrType = event[isrType].mother1();
              if (isrType > 2) isrType -= beamOffset;
              dipNow.isrType   = isrType;
              done             = true;
            }
          } 
          // If above options failed, then create new dipole
          if (!done) {
            int iRadNow = partonSystemsPtr->getIndexOfOut(dipNow.system, iRad);
            if (iRadNow != -1)
              setupQEDdip(dipNow.system, iRadNow, dipNow.chgType,
                          dipNow.gamType, event, true);
            else
              infoPtr->errorMsg("Warning in TimeShower::rescatterUpdate: "
              "failed to locate radiator in system");

            dipNow.colType = 0;
            dipNow.chgType = 0;
            dipNow.gamType = 0;

            infoPtr->errorMsg("Warning in TimeShower::rescatterUpdate: "
            "failed to locate new recoiling charge partner");
          }
        }
      }

    // End of loop over dipoles and two incoming sides. 
    }
  }
   
}

//--------------------------------------------------------------------------

// Update dipole list after each ISR emission (so not used for resonances).  

void TimeShower::update( int iSys, Event& event) {

  // Start list of rescatterers that gave further changed systems in ISR.
  vector<int> iRescatterer;

  // Find new and old positions of incoming partons in the system.
  vector<int> iNew, iOld;
  iNew.push_back( partonSystemsPtr->getInA(iSys) );
  iOld.push_back( event[iNew[0]].daughter2() ); 
  iNew.push_back( partonSystemsPtr->getInB(iSys) );
  iOld.push_back( event[iNew[1]].daughter2() ); 

  // Ditto for outgoing partons, except the newly created one.
  int sizeOut = partonSystemsPtr->sizeOut(iSys) - 1;  
  for (int i = 0; i < sizeOut; ++i) {
    int iNow = partonSystemsPtr->getOut(iSys, i);
    iNew.push_back( iNow );
    iOld.push_back( event[iNow].mother1() );
    // Add non-final to list of rescatterers.
    if (!event[iNow].isFinal()) iRescatterer.push_back( iNow );
  }
  int iNewNew = partonSystemsPtr->getOut(iSys, sizeOut);
  
  // Swap beams to let 0 be side on which branching occured. 
  if (event[iNew[0]].status() != -41) {
    swap( iNew[0], iNew[1]);   
    swap( iOld[0], iOld[1]);   
  }

  // Loop over all dipole ends belonging to the system 
  // or to the recoil system, if different.
  for (int iDip = 0; iDip < int(dipEnd.size()); ++iDip) 
  if (dipEnd[iDip].system == iSys || dipEnd[iDip].systemRec == iSys) {
    TimeDipoleEnd& dipNow = dipEnd[iDip];

    // Replace radiator (always in final state so simple).
    for (int i = 2; i < 2 + sizeOut; ++i) 
    if (dipNow.iRadiator == iOld[i]) { 
      dipNow.iRadiator = iNew[i];
      break;
    }

    // Replace ME partner (always in final state, if exists, so simple).
    for (int i = 2; i < 2 + sizeOut; ++i) 
    if (dipNow.iMEpartner == iOld[i]) { 
      dipNow.iMEpartner = iNew[i];
      break;
    }

    // Recoiler: by default pick old one, only moved. Note excluded beam.
    int iRec = 0;
    if (dipNow.systemRec == iSys) { 
      for (int i = 1; i < 2 + sizeOut; ++i) 
      if (dipNow.iRecoiler == iOld[i]) { 
        iRec = iNew[i];
        break;
      }

      // QCD recoiler: check if colour hooks up with new final parton.
      if ( dipNow.colType > 0 
        && event[dipNow.iRadiator].col() == event[iNewNew].acol() ) { 
        iRec = iNewNew;
        dipNow.isrType = 0;
      }
      if ( dipNow.colType < 0 
        && event[dipNow.iRadiator].acol() == event[iNewNew].col() ) { 
        iRec = iNewNew;
        dipNow.isrType = 0;
      }
      
      // QCD recoiler: check if colour hooks up with new beam parton.
      if ( iRec == 0 && dipNow.colType > 0 
        && event[dipNow.iRadiator].col()  == event[iNew[0]].col() ) 
        iRec = iNew[0];
      if ( iRec == 0 && dipNow.colType < 0 
        && event[dipNow.iRadiator].acol() == event[iNew[0]].acol() )   
        iRec = iNew[0];

      // QED/photon recoiler: either to new particle or remains to beam.
      if ( iRec == 0 && (dipNow.chgType != 0 || dipNow.gamType != 0) ) {
        if ( event[iNew[0]].chargeType() == 0 ) {
          iRec = iNewNew;
          dipNow.isrType = 0;
        } else {
          iRec = iNew[0];
        }
      }

    // Recoiler in another system: keep it as is.
    } else iRec = dipNow.iRecoiler;

    // Done. Kill dipole if failed to find new recoiler. 
    dipNow.iRecoiler = iRec;
    if ( iRec == 0 && (dipNow.colType != 0 || dipNow.chgType != 0
      || dipNow.gamType != 0) ) {
      dipNow.colType = 0;
      dipNow.chgType = 0;
      dipNow.gamType = 0; 
      infoPtr->errorMsg("Error in TimeShower::update: "
      "failed to locate new recoiling partner");
    } 
  }
  
  // Find new dipole end formed by colour index.
  int colTag = event[iNewNew].col();     
  if (doQCDshower && colTag > 0) 
    setupQCDdip( iSys, sizeOut, colTag, 1, event, false, true); 
  
  // Find new dipole end formed by anticolour index.
  int acolTag = event[iNewNew].acol();     
  if (doQCDshower && acolTag > 0) 
    setupQCDdip( iSys, sizeOut, acolTag, -1, event, false, true); 

  // Find new "charge-dipole" and "photon-dipole" ends. 
  int  chgType  = event[iNewNew].chargeType();  
  bool doChgDip = (chgType != 0) 
                  && ( ( doQEDshowerByQ && event[iNewNew].isQuark()  )
                    || ( doQEDshowerByL && event[iNewNew].isLepton() ) );
  int  gamType  = (event[iNewNew].id() == 22) ? 1 : 0;
  bool doGamDip = (gamType == 1) && doQEDshowerByGamma;
  if (doChgDip || doGamDip) 
    setupQEDdip( iSys, sizeOut, chgType, gamType, event, true); 

  // Start iterate over list of rescatterers - may be empty.
  int iRescNow = -1;
  while (++iRescNow < int(iRescatterer.size())) {

    // Identify systems that rescatterers belong to.
    int iOutNew    = iRescatterer[iRescNow];    
    int iInNew     = event[iOutNew].daughter1();
    int iSysResc   = partonSystemsPtr->getSystemOf(iInNew, true);

    // Find new and old positions of incoming partons in the system.
    iNew.resize(0); 
    iOld.resize(0);
    iNew.push_back( partonSystemsPtr->getInA(iSysResc) );
    iOld.push_back( event[iNew[0]].daughter1() ); 
    iNew.push_back( partonSystemsPtr->getInB(iSysResc) );
    iOld.push_back( event[iNew[1]].daughter1() ); 

    // Ditto for outgoing partons.
    sizeOut = partonSystemsPtr->sizeOut(iSysResc);  
    for (int i = 0; i < sizeOut; ++i) {
      int iNow = partonSystemsPtr->getOut(iSysResc, i);
      iNew.push_back( iNow );
      iOld.push_back( event[iNow].mother1() );
      // Add non-final to list of rescatterers.
      if (!event[iNow].isFinal()) iRescatterer.push_back( iNow );
    }

    // Loop over all dipole ends belonging to the system 
    // or to the recoil system, if different.
    for (int iDip = 0; iDip < int(dipEnd.size()); ++iDip) 
    if (dipEnd[iDip].system == iSysResc 
      || dipEnd[iDip].systemRec == iSysResc) {
      TimeDipoleEnd& dipNow = dipEnd[iDip];

      // Replace radiator (always in final state so simple).
      for (int i = 2; i < 2 + sizeOut; ++i) 
      if (dipNow.iRadiator == iOld[i]) { 
        dipNow.iRadiator = iNew[i];
        break;
      }

      // Replace ME partner (always in final state, if exists, so simple).
      for (int i = 2; i < 2 + sizeOut; ++i) 
      if (dipNow.iMEpartner == iOld[i]) { 
        dipNow.iMEpartner = iNew[i];
        break;
      }

      // Replace recoiler.
      for (int i = 0; i < 2 + sizeOut; ++i) 
      if (dipNow.iRecoiler == iOld[i]) { 
        dipNow.iRecoiler = iNew[i];
        break;
      }
    }

  // End iterate over list of rescatterers.
  }

}

//--------------------------------------------------------------------------

// Setup a dipole end for a QCD colour charge.

void TimeShower::setupQCDdip( int iSys, int i, int colTag, int colSign, 
  Event& event, bool isOctetOnium, bool limitPTmaxIn) {
 
  // Initial values. Find if allowed to hook up beams.
  int iRad     = partonSystemsPtr->getOut(iSys, i);
  int iRec     = 0;
  int sizeAllA = partonSystemsPtr->sizeAll(iSys);
  int sizeOut  = partonSystemsPtr->sizeOut(iSys);
  int sizeAll  = ( allowBeamRecoil ) ? sizeAllA : sizeOut;
  int sizeIn   = sizeAll - sizeOut;
  int sizeInA  = sizeAllA - sizeIn - sizeOut;
  int iOffset  = i + sizeAllA - sizeOut;
  bool otherSystemRec = false;
  bool allowInitial   = (partonSystemsPtr->hasInAB(iSys)) ? true : false;
  // PS dec 2010: possibility to allow for several recoilers and each with
  // flexible normalization
  bool   isFlexible   = false;
  double flexFactor   = 1.0;
  vector<int> iRecVec(0);

  // Colour: other end by same index in beam or opposite in final state.
  // Exclude rescattered incoming and not final outgoing.
  if (colSign > 0) 
  for (int j = 0; j < sizeAll; ++j) if (j + sizeInA != iOffset) {
    int iRecNow = partonSystemsPtr->getAll(iSys, j + sizeInA);
    if ( ( j <  sizeIn && event[iRecNow].col()  == colTag
      && !event[iRecNow].isRescatteredIncoming() )
      || ( j >= sizeIn && event[iRecNow].acol() == colTag 
      && event[iRecNow].isFinal() ) ) { 
      iRec = iRecNow;
      break;
    }
  }
 
  // Anticolour: other end by same index in beam or opposite in final state.
  // Exclude rescattered incoming and not final outgoing.
  if (colSign < 0) 
  for (int j = 0; j < sizeAll; ++j) if (j + sizeInA != iOffset) {
    int iRecNow = partonSystemsPtr->getAll(iSys, j + sizeInA);
    if ( ( j <  sizeIn && event[iRecNow].acol()  == colTag
      && !event[iRecNow].isRescatteredIncoming() )
      || ( j >= sizeIn && event[iRecNow].col() == colTag
      && event[iRecNow].isFinal() ) ) { 
      iRec = iRecNow;
      break;
    }
  }

  // Resonance decays (= no instate): 
  // other end to nearest recoiler in same system final state,
  // by (p_i + p_j)^2 - (m_i + m_j)^2 = 2 (p_i p_j - m_i m_j).  
  // (junction colours more involved, so keep track if junction colour)
  bool hasJunction = false;
  if (iRec == 0 && !allowInitial) {
    for (int iJun = 0; iJun < event.sizeJunction(); ++ iJun) {
      // For types 1&2, all legs in final state
      // For types 3&4, two legs in final state
      // For types 5&6, one leg in final state
      int iBeg = (event.kindJunction(iJun)-1)/2;
      for (int iLeg = iBeg; iLeg < 3; ++iLeg) 
	if (event.endColJunction( iJun, iLeg) == colTag) hasJunction  = true; 
    }
    double ppMin = LARGEM2; 
    for (int j = 0; j < sizeOut; ++j) if (j != i) { 
        int iRecNow  = partonSystemsPtr->getOut(iSys, j);
        if (!event[iRecNow].isFinal()) continue;
        double ppNow = event[iRecNow].p() * event[iRad].p() 
	  - event[iRecNow].m() * event[iRad].m();
        if (ppNow < ppMin) {
          iRec  = iRecNow;
          ppMin = ppNow;
        } 
      }
  }

  // If no success then look for matching (anti)colour anywhere in final state.
  if ( iRec == 0 || (!doInterleave && !event[iRec].isFinal()) ) {
    iRec = 0;
    for (int j = 0; j < event.size(); ++j) if (event[j].isFinal()) 
    if ( (colSign > 0 && event[j].acol() == colTag)
      || (colSign < 0 && event[j].col()  == colTag) ) {
      iRec = j;
      otherSystemRec = true;
      break;
    }

    // If no success then look for match to non-rescattered in initial state.
    if (iRec == 0 && allowInitial) {
      for (int iSysR = 0; iSysR < partonSystemsPtr->sizeSys(); ++iSysR) 
      if (iSysR != iSys) {
        int j = partonSystemsPtr->getInA(iSysR);
        if (j > 0 && event[j].isRescatteredIncoming()) j = 0;
        if (j > 0 && ( (colSign > 0 && event[j].col() == colTag)
          || (colSign < 0 && event[j].acol()  == colTag) ) ) {
          iRec = j;
          otherSystemRec = true;
          break;
        }
        j = partonSystemsPtr->getInB(iSysR);
        if (j > 0 && event[j].isRescatteredIncoming()) j = 0;
        if (j > 0 && ( (colSign > 0 && event[j].col() == colTag)
          || (colSign < 0 && event[j].acol()  == colTag) ) ) {
          iRec = j;
          otherSystemRec = true;
          break;
        }
      }
    }
  }

  // Junctions (PS&ND dec 2010)
  // For types 1&2: all legs in final state
  //                half-strength dipoles between all legs
  // For types 3&4, two legs in final state
  //                full-strength dipole between final-state legs 
  // For types 5&6, one leg in final state	
  //                no final-state dipole end
  
  if (hasJunction) {
    for (int iJun = 0; iJun < event.sizeJunction(); ++ iJun) {
      int kindJun = event.kindJunction(iJun);
      int iBeg = (kindJun-1)/2;
      for (int iLeg = iBeg; iLeg < 3; ++iLeg) {	
	if (event.endColJunction( iJun, iLeg) == colTag) {
	  // For types 5&6, no other leg to recoil against. Switch off if
	  // no other particles at all, since radiation then handled by ISR. 
	  // Example: qq -> ~t* : no radiation off ~t*
	  // Allow radiation + recoil if unconnected partners available
	  // Example: qq -> ~t* -> tbar ~chi0 : allow radiation off tbar, 
	  //                                    with ~chi0 as recoiler
	  if (kindJun >= 5) { 
	    if (sizeOut == 1) return;
	    else break;
	  }
	  // For junction types 3 & 4, span one full-strength dipole
	  // (only look inside same decay system)
	  else if (kindJun >= 3) {
	    int iLegRec = 3-iLeg;
	    int colTagRec = event.endColJunction( iJun, iLegRec);
	    for (int j = 0; j < sizeOut; ++j) if (j != i) { 
		int iRecNow  = partonSystemsPtr->getOut(iSys, j);
		if (!event[iRecNow].isFinal()) continue;
		if ( (colSign > 0 && event[iRecNow].col()  == colTagRec)
		  || (colSign < 0 && event[iRecNow].acol() == colTagRec) ) {
		  // Only accept if staying inside same system
		  iRec = iRecNow;		  
		  break;
		}
	      }
	  }
	  // For junction types 1 & 2, span two half-strength dipoles
	  // (only look inside same decay system)
	  else {
	    // Loop over two half-strength dipole connections
	    for (int jLeg = 1; jLeg <= 2; jLeg++) {
	      int iLegRec = (iLeg + jLeg) % 3;
	      int colTagRec = event.endColJunction( iJun, iLegRec);
	      for (int j = 0; j < sizeOut; ++j) if (j != i) { 
		  int iRecNow  = partonSystemsPtr->getOut(iSys, j);
		  if (!event[iRecNow].isFinal()) continue;
		  if ( (colSign > 0 && event[iRecNow].col()  == colTagRec)
		    || (colSign < 0 && event[iRecNow].acol() == colTagRec) ) {
		    // Store recoilers in temporary array
		    iRecVec.push_back(iRecNow);
		    // Set iRec != 0 for checks below
		    iRec = iRecNow;
		  }
		}
	    }
	    
	  }     // End if-then-else of junction kinds
	  
	}       // End if leg has right color tag 
      }         // End of loop over junction legs
    }           // End loop over junctions
    
  }             // End main junction if
  
  // If fail, then other end to nearest recoiler in same system final state,
  // by (p_i + p_j)^2 - (m_i + m_j)^2 = 2 (p_i p_j - m_i m_j).
  if (iRec == 0) {    
    double ppMin = LARGEM2; 
    for (int j = 0; j < sizeOut; ++j) if (j != i) { 
      int iRecNow  = partonSystemsPtr->getOut(iSys, j);
      if (!event[iRecNow].isFinal()) continue;
      double ppNow = event[iRecNow].p() * event[iRad].p() 
                   - event[iRecNow].m() * event[iRad].m();
      if (ppNow < ppMin) {
        iRec  = iRecNow;
        ppMin = ppNow;
      } 
    }     
  }  

  // If fail, then other end to nearest recoiler in any system final state,
  // by (p_i + p_j)^2 - (m_i + m_j)^2 = 2 (p_i p_j - m_i m_j).
  if (iRec == 0) {
    double ppMin = LARGEM2; 
    for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow) 
    if (iRecNow != iRad && event[iRecNow].isFinal()) {
      double ppNow = event[iRecNow].p() * event[iRad].p() 
                   - event[iRecNow].m() * event[iRad].m();
      if (ppNow < ppMin) {
        iRec  = iRecNow;
        otherSystemRec = true;
        ppMin = ppNow;
      } 
    }     
  }  

  // PS dec 2010: make sure iRec is stored in iRecVec 
  if (iRecVec.size() == 0 && iRec != 0) iRecVec.push_back(iRec);
    
  // Remove any zero recoilers from normalization
  int nRec = iRecVec.size();
  for (unsigned int mRec = 0; mRec < iRecVec.size(); ++mRec) 
    if (iRecVec[mRec] <= 0) nRec--;
  if (nRec >= 2) {
    isFlexible = true;
    flexFactor = 1.0/nRec;
  }
  
  // Check for failure to locate any recoiler
  if ( nRec <= 0 ) { 
    infoPtr->errorMsg("Error in TimeShower::setupQCDdip: "
		      "failed to locate any recoiling partner");
    return;
  }
  
  // Store dipole colour end(s).
  for (unsigned int mRec = 0; mRec < iRecVec.size(); ++mRec) {
    iRec = iRecVec[mRec];
    if (iRec <= 0) continue;
    // Max scale either by parton scale or by half dipole mass.
    double pTmax = event[iRad].scale();
    if (limitPTmaxIn) {
      if (iSys == 0 || (iSys == 1 && doSecondHard)) pTmax *= pTmaxFudge;
      else if (sizeIn > 0) pTmax *= pTmaxFudgeMPI;
    } else pTmax = 0.5 * m( event[iRad], event[iRec]);
    int colType  = (event[iRad].id() == 21) ? 2 * colSign : colSign;
    int isrType  = (event[iRec].isFinal()) ? 0 : event[iRec].mother1();
    // This line in case mother is a rescattered parton.
    while (isrType > 2 + beamOffset) isrType = event[isrType].mother1();
    if (isrType > 2) isrType -= beamOffset;
    dipEnd.push_back( TimeDipoleEnd( iRad, iRec, pTmax, 
      colType, 0, 0, isrType, iSys, -1, -1, isOctetOnium) );

    // If hooked up with other system then find which.
    if (otherSystemRec) {
      int systemRec = partonSystemsPtr->getSystemOf(iRec, true);
      if (systemRec >= 0) dipEnd.back().systemRec = systemRec;
      dipEnd.back().MEtype = 0;
    } 

    // PS dec 2010
    // If non-unity (flexible) normalization, set normalization factor
    if (isFlexible) {
      dipEnd.back().isFlexible = true;
      dipEnd.back().flexFactor = flexFactor;
    }
  } 

}

//--------------------------------------------------------------------------

// Setup a dipole end for a QED colour charge or a photon.
// No failsafe choice of recoiler, so gradually widen search.

void TimeShower::setupQEDdip( int iSys, int i, int chgType, int gamType,
  Event& event, bool limitPTmaxIn) {

  // Initial values. Find if allowed to hook up beams.
  int iRad     = partonSystemsPtr->getOut(iSys, i);
  int idRad    = event[iRad].id();
  int iRec     = 0;
  int sizeAllA = partonSystemsPtr->sizeAll(iSys);
  int sizeOut  = partonSystemsPtr->sizeOut(iSys);
  int sizeAll  = ( allowBeamRecoil ) ? sizeAllA : sizeOut;
  int sizeIn   = sizeAll - sizeOut;
  int sizeInA  = sizeAllA - sizeIn - sizeOut;
  int iOffset  = i + sizeAllA - sizeOut;
  double ppMin = LARGEM2;
  bool hasRescattered = false;
  bool otherSystemRec = false;

  // Find nearest same- (opposide-) flavour recoiler in initial (final)
  // state of same system, excluding rescattered (in or out) partons.
  // Also find if system is involved in rescattering.
  // Note: (p_i + p_j)2 - (m_i + m_j)2 = 2 (p_i p_j - m_i m_j).
  for (int j = 0; j < sizeAll; ++j) if (j + sizeInA != iOffset) {
    int iRecNow  = partonSystemsPtr->getAll(iSys, j + sizeInA);
    if ( (j <  sizeIn && !event[iRecNow].isRescatteredIncoming())
      || (j >= sizeIn && event[iRecNow].isFinal()) ) {
      if ( (j <  sizeIn && event[iRecNow].id() ==  idRad)
        || (j >= sizeIn && event[iRecNow].id() == -idRad) ) {
        double ppNow = event[iRecNow].p() * event[iRad].p()
                     - event[iRecNow].m() * event[iRad].m();
        if (ppNow < ppMin) {
          iRec  = iRecNow;
          ppMin = ppNow;
        }
      }
    } else hasRescattered = true;
  }

  // If rescattering then find nearest opposite-flavour recoiler 
  // anywhere in final state.
  if (iRec == 0 && hasRescattered) {
    for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow)
    if (event[iRecNow].id() == -idRad && event[iRecNow].isFinal()) {
      double ppNow = event[iRecNow].p() * event[iRad].p()
                   - event[iRecNow].m() * event[iRad].m();
      if (ppNow < ppMin) {
        iRec  = iRecNow;
        ppMin = ppNow;
        otherSystemRec = true;
      }
    }
  }

  // Find nearest recoiler in same system, charge-squared-weighted,
  // including initial state, but excluding rescatterer.
  if (iRec == 0)
  for (int j = 0; j < sizeAll; ++j) if (j + sizeInA != iOffset) {
    int iRecNow = partonSystemsPtr->getAll(iSys, j + sizeInA);
    int chgTypeRecNow = event[iRecNow].chargeType();
    if (chgTypeRecNow == 0) continue;
    if ( (j <  sizeIn && !event[iRecNow].isRescatteredIncoming())
      || (j >= sizeIn && event[iRecNow].isFinal()) ) {
      double ppNow = (event[iRecNow].p() * event[iRad].p()
                   -  event[iRecNow].m() * event[iRad].m())
                   / pow2(chgTypeRecNow);
      if (ppNow < ppMin) {
        iRec  = iRecNow;
        ppMin = ppNow;
      }
    }
  }

  // If rescattering then find nearest recoiler in the final state, 
  // charge-squared-weighted.
  if (iRec == 0 && hasRescattered) {
    for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow)
    if (iRecNow != iRad && event[iRecNow].isFinal()) {
      int chgTypeRecNow = event[iRecNow].chargeType();
      if (chgTypeRecNow != 0 && event[iRecNow].isFinal()) {
        double ppNow = (event[iRecNow].p() * event[iRad].p()
                     -  event[iRecNow].m() * event[iRad].m())
                     / pow2(chgTypeRecNow);
        if (ppNow < ppMin) {
          iRec  = iRecNow;
          ppMin = ppNow;
          otherSystemRec = true;
        }
      }
    }
  }

  // Find any nearest recoiler in final state of same system.
  if (iRec == 0)
  for (int j = 0; j < sizeOut; ++j) if (j != i) {
    int iRecNow  = partonSystemsPtr->getOut(iSys, j);
    double ppNow = event[iRecNow].p() * event[iRad].p()
                 - event[iRecNow].m() * event[iRad].m();
    if (ppNow < ppMin) {
      iRec  = iRecNow;
      ppMin = ppNow;
    }
  }

  // Find any nearest recoiler in final state.
  if (iRec == 0)
  for (int iRecNow = 0; iRecNow < event.size(); ++iRecNow)
  if (iRecNow != iRad && event[iRecNow].isFinal()) {
    double ppNow = event[iRecNow].p() * event[iRad].p()
                 - event[iRecNow].m() * event[iRad].m();
    if (ppNow < ppMin) {
      iRec  = iRecNow;
      ppMin = ppNow;
      otherSystemRec = true;
    }
  }

  // Fill charge-dipole or photon-dipole end.
  if (iRec > 0) {
    // Max scale either by parton scale or by half dipole mass.
    double pTmax = event[iRad].scale();
    if (limitPTmaxIn) {
      if (iSys == 0 || (iSys == 1 && doSecondHard)) pTmax *= pTmaxFudge;
      else if (sizeIn > 0) pTmax *= pTmaxFudgeMPI;
    } else pTmax = 0.5 * m( event[iRad], event[iRec]);
    int isrType = (event[iRec].isFinal()) ? 0 : event[iRec].mother1();
    // This line in case mother is a rescattered parton.
    while (isrType > 2 + beamOffset) isrType = event[isrType].mother1();
    if (isrType > 2) isrType -= beamOffset;
    dipEnd.push_back( TimeDipoleEnd(iRad, iRec, pTmax,
      0, chgType, gamType, isrType, iSys, -1) );

    // If hooked up with other system then find which.
    if (otherSystemRec) {
      int systemRec = partonSystemsPtr->getSystemOf(iRec);
      if (systemRec >= 0) dipEnd.back().systemRec = systemRec;
      dipEnd.back().MEtype = 0;
    }

  // Failure to find other end of dipole.
  } else {
    infoPtr->errorMsg("Error in TimeShower::setupQEDdip: "
      "failed to locate any recoiling partner");
  }

}

//--------------------------------------------------------------------------

// Setup a dipole end for a Hidden Valley colour charge.

void TimeShower::setupHVdip( int iSys, int i, Event& event, 
  bool limitPTmaxIn) {
 
  // Initial values.
  int iRad    = partonSystemsPtr->getOut(iSys, i);
  int iRec    = 0;
  int idRad   = event[iRad].id();
  int sizeOut = partonSystemsPtr->sizeOut(iSys);

  // Hidden Valley colour positive for positive id, and vice versa.
  // Find opposte HV colour in final state of same system.
  for (int j = 0; j < sizeOut; ++j) if (j != i) {
    int iRecNow = partonSystemsPtr->getOut(iSys, j);
    int idRec   = event[iRecNow].id();
    if ( (abs(idRec) > 4900000 && abs(idRec) < 4900017)
      && idRad * idRec < 0) {
      iRec = iRecNow;
      break;
    }
  }

  // Else find heaviest other final-state in same system.
  // (Intended for decays; should mainly be two-body so unique.)
  double mMax = -sqrt(LARGEM2);
   if (iRec == 0)
  for (int j = 0; j < sizeOut; ++j) if (j != i) {
    int iRecNow = partonSystemsPtr->getOut(iSys, j);
    if (event[iRecNow].m() > mMax) {
      iRec = iRecNow;
      mMax = event[iRecNow].m();
    }
  }

  // Set up dipole end, or report failure. 
  if (iRec > 0) {
    // Max scale either by parton scale or by half dipole mass.
    double pTmax = event[iRad].scale();
    if (limitPTmaxIn) {
      if (iSys == 0 || (iSys == 1 && doSecondHard)) pTmax *= pTmaxFudge;
    } else pTmax = 0.5 * m( event[iRad], event[iRec]);
    int colvType  = (event[iRad].id() > 0) ? 1 : -1;
    dipEnd.push_back( TimeDipoleEnd( iRad, iRec, pTmax, 0, 0, 0, 0, 
      iSys, -1, -1, false, true, colvType) );
  } else infoPtr->errorMsg("Error in TimeShower::setupHVdip: "
      "failed to locate any recoiling partner");

}
 
//--------------------------------------------------------------------------

// Select next pT in downwards evolution of the existing dipoles.

double TimeShower::pTnext( Event& event, double pTbegAll, double pTendAll) {

  // Begin loop over all possible radiating dipole ends.
  dipSel  = 0;
  iDipSel = -1;
  double pT2sel = pTendAll * pTendAll;
  for (int iDip = 0; iDip < int(dipEnd.size()); ++iDip) {
    TimeDipoleEnd& dip = dipEnd[iDip]; 

    // Check if global recoil should be used.
    useLocalRecoilNow = !(globalRecoil && dip.system == 0 
      && partonSystemsPtr->sizeOut(0) <= nMaxGlobalRecoil);
   
    // Dipole properties; normal local recoil. 
    dip.mRad   = event[dip.iRadiator].m(); 
    if (useLocalRecoilNow) {
      dip.mRec = event[dip.iRecoiler].m(); 
      dip.mDip = m( event[dip.iRadiator], event[dip.iRecoiler] );

    // Dipole properties, alternative global recoil. Squares.
    } else {
      Vec4 pSumGlobal;
      for (int i = 0; i < partonSystemsPtr->sizeOut( dip.system); ++i) { 
        int ii = partonSystemsPtr->getOut( dip.system, i);
        if (ii !=  dip.iRadiator) pSumGlobal += event[ii].p();
      }
      dip.mRec = pSumGlobal.mCalc();
      dip.mDip = m( event[dip.iRadiator].p(), pSumGlobal); 
    } 
    dip.m2Rad  = pow2(dip.mRad);
    dip.m2Rec  = pow2(dip.mRec);
    dip.m2Dip  = pow2(dip.mDip);

    // Find maximum evolution scale for dipole.
    dip.m2DipCorr    = pow2(dip.mDip - dip.mRec) - dip.m2Rad; 
    double pTbegDip  = min( pTbegAll, dip.pTmax ); 
    double pT2begDip = min( pow2(pTbegDip), 0.25 * dip.m2DipCorr);

    // Do QCD, QED or HV evolution if it makes sense.
    dip.pT2 = 0.;
    if (pT2begDip > pT2sel) {
      if      (dip.colType != 0) 
        pT2nextQCD(pT2begDip, pT2sel, dip, event);
      else if (dip.chgType != 0 || dip.gamType != 0)                 
        pT2nextQED(pT2begDip, pT2sel, dip, event);
      else if (dip.colvType != 0)
        pT2nextHV(pT2begDip, pT2sel, dip, event);

      // Update if found larger pT than current maximum.
      if (dip.pT2 > pT2sel) {
        pT2sel  = dip.pT2;
        dipSel  = &dip;
        iDipSel = iDip;
      }
    } 
  } 

  // Return nonvanishing value if found pT bigger than already found.
  return (dipSel == 0) ? 0. : sqrt(pT2sel); 

}

//--------------------------------------------------------------------------

// Evolve a QCD dipole end. 

void TimeShower::pT2nextQCD(double pT2begDip, double pT2sel, 
  TimeDipoleEnd& dip, Event& event) { 

  // Lower cut for evolution. Return if no evolution range.
  double pT2endDip = max( pT2sel, pT2colCut );   
  if (pT2begDip < pT2endDip) return;   

  // Upper estimate for matrix element weighting and colour factor.
  // Special cases for triplet recoiling against gluino and octet onia.
  // Note that g -> g g and g -> q qbar are split on two sides.
  int    colTypeAbs = abs(dip.colType);
  double wtPSglue   = 2.;
  double colFac     = (colTypeAbs == 1) ? 4./3. : 3./2.;
  if (dip.MEgluinoRec)  colFac  = 3.;  
  if (dip.isOctetOnium) colFac *= 0.5 * octetOniumColFac;
  // PS dec 2010. Include possibility for flexible normalization,
  // e.g., for dipoles stretched to junctions or to switch off radiation.
  if (dip.isFlexible)   colFac *= dip.flexFactor;
  double wtPSqqbar  = (colTypeAbs == 2) ? 0.25 * nGluonToQuark : 0.;
  
  // Variables used inside evolution loop. (Mainly dummy start values.)
  dip.pT2              = pT2begDip;
  int    nFlavour      = 3;
  double zMinAbs       = 0.5;
  double pT2min        = pT2endDip;
  double b0            = 4.5;
  double Lambda2       = Lambda3flav2; 
  double emitCoefGlue  = 0.;
  double emitCoefQqbar = 0.; 
  double emitCoefTot   = 0.; 
  double wt            = 0.; 
  bool   mustFindRange = true;
  
  // Begin evolution loop towards smaller pT values.
  do { 

    // Initialize evolution coefficients at the beginning and
    // reinitialize when crossing c and b flavour thresholds.
    if (mustFindRange) {

      // Determine overestimated z range; switch at c and b masses.
      if (dip.pT2 > m2b) {
        nFlavour = 5;
        pT2min   = max( m2b, pT2endDip); 
        b0       = 23./6.;
        Lambda2  = Lambda5flav2;
      } else if (dip.pT2 > m2c) {
        nFlavour = 4;
        pT2min   = max( m2c, pT2endDip); 
        b0       = 25./6.;
        Lambda2  = Lambda4flav2;
      } else { 
        nFlavour = 3;
        pT2min   = pT2endDip;
        b0       = 27./6.;
        Lambda2  = Lambda3flav2;
      }
      // A change of renormalization scale expressed by a change of Lambda. 
      Lambda2 /= renormMultFac;
      zMinAbs = 0.5 - sqrtpos( 0.25 - pT2min / dip.m2DipCorr );
      if (zMinAbs < SIMPLIFYROOT) zMinAbs = pT2min / dip.m2DipCorr;

      // Find emission coefficients for X -> X g and g -> q qbar.
      emitCoefGlue = wtPSglue * colFac * log(1. / zMinAbs - 1.);
      emitCoefTot  = emitCoefGlue;
      if (colTypeAbs == 2 && event[dip.iRadiator].id() == 21) {
        emitCoefQqbar = wtPSqqbar * (1. - 2. * zMinAbs);
        emitCoefTot  += emitCoefQqbar;
      }

      // Initialization done for current range.
      mustFindRange = false;
    } 

    // Pick pT2 (in overestimated z range) for fixed alpha_strong.
    if (alphaSorder == 0) {
      dip.pT2 = dip.pT2 * pow( rndmPtr->flat(), 
        1. / (alphaS2pi * emitCoefTot) );

    // Ditto for first-order alpha_strong.
    } else if (alphaSorder == 1) {
      dip.pT2 = Lambda2 * pow( dip.pT2 / Lambda2, 
        pow( rndmPtr->flat(), b0 / emitCoefTot) );

      // For second order reject by second term in alpha_strong expression.
    } else {
      do dip.pT2 = Lambda2 * pow( dip.pT2 / Lambda2, 
        pow( rndmPtr->flat(), b0 / emitCoefTot) );
      while (alphaS.alphaS2OrdCorr(renormMultFac * dip.pT2) < rndmPtr->flat() 
        && dip.pT2 > pT2min);
    }
    wt = 0.;
  
    // If crossed c or b thresholds: continue evolution from threshold.
    if (nFlavour == 5 && dip.pT2 < m2b) {  
      mustFindRange = true;
      dip.pT2       = m2b;
    } else if ( nFlavour == 4 && dip.pT2 < m2c) { 
      mustFindRange = true;
      dip.pT2       = m2c;

    // Abort evolution if below cutoff scale, or below another branching.
    } else {
      if ( dip.pT2 < pT2endDip) { dip.pT2 = 0.; return; }

      // Pick kind of branching: X -> X g or g -> q qbar. 
      dip.flavour  = 21;
      dip.mFlavour = 0.;
      if (colTypeAbs == 2 && emitCoefQqbar > rndmPtr->flat() 
        * emitCoefTot) dip.flavour = 0; 

      // Pick z: either dz/(1-z) or flat dz.
      if (dip.flavour == 21) {
        dip.z = 1. - zMinAbs * pow( 1. / zMinAbs - 1., rndmPtr->flat() );
      } else { 
        dip.z = zMinAbs + (1. - 2. * zMinAbs) * rndmPtr->flat();   
      }
  
      // Do not accept branching if outside allowed z range.
      double zMin = 0.5 - sqrtpos( 0.25 - dip.pT2 / dip.m2DipCorr ); 
      if (zMin < SIMPLIFYROOT) zMin = dip.pT2 / dip.m2DipCorr;
      dip.m2 = dip.m2Rad + dip.pT2 / (dip.z * (1. - dip.z));
      if (dip.z > zMin && dip.z < 1. - zMin 
        && dip.m2 * dip.m2Dip < dip.z * (1. - dip.z) 
        * pow2(dip.m2Dip + dip.m2 - dip.m2Rec) ) {

        // Flavour choice for g -> q qbar.
        if (dip.flavour == 0) {
          dip.flavour  = min(5, 1 + int(nGluonToQuark * rndmPtr->flat())); 
          dip.mFlavour = particleDataPtr->m0(dip.flavour);
        }

        // No z weight, except threshold, if to do ME corrections later on.
        if (dip.MEtype > 0) { 
          wt = 1.;
          if (dip.flavour < 10 && dip.m2 < THRESHM2 * pow2(dip.mFlavour))
            wt = 0.; 

        // z weight for X -> X g.
        } else if (dip.flavour == 21 
          && (colTypeAbs == 1 || colTypeAbs == 3) ) {
          wt = (1. + pow2(dip.z)) / wtPSglue;
        } else if (dip.flavour == 21) {     
          wt = (1. + pow3(dip.z)) / wtPSglue;
           
        // z weight for g -> q qbar.
        } else {
          double beta  = sqrtpos( 1. - 4. * pow2(dip.mFlavour) / dip.m2 );
          wt = beta * ( pow2(dip.z) + pow2(1. - dip.z) );
        }

        // Suppression factors for dipole to beam remnant.
        if (dip.isrType != 0 && useLocalRecoilNow) {
         BeamParticle& beam = (dip.isrType == 1) ? *beamAPtr : *beamBPtr;
          int iSysRec = dip.systemRec;
          double xOld = beam[iSysRec].x();
          double xNew = xOld * (1. + (dip.m2 - dip.m2Rad) / 
            (dip.m2Dip - dip.m2Rad));
          double xMaxAbs = beam.xMax(iSysRec);
          if (xMaxAbs < 0.) {
            infoPtr->errorMsg("Warning in TimeShower::pT2nextQCD: "
            "xMaxAbs negative"); 
            return;
          }
 
          // Firstly reduce by PDF ratio.
          if (xNew > xMaxAbs) wt = 0.;              
          else {
            int idRec     = event[dip.iRecoiler].id();
            double pdfOld = max ( TINYPDF, 
              beam.xfISR( iSysRec, idRec, xOld, factorMultFac * dip.pT2) ); 
            double pdfNew = 
              beam.xfISR( iSysRec, idRec, xNew, factorMultFac * dip.pT2); 
            wt *= min( 1., pdfNew / pdfOld); 
          }

          // Secondly optionally reduce by 4 pT2_hard / (4 pT2_hard + m2).
          if (dampenBeamRecoil) {
            double pTpT = sqrt(event[dip.iRadiator].pT2() * dip.pT2);
            wt *= pTpT / (pTpT + dip.m2);
          }
        }

        // Optional dampening of large pT values in hard system.
        if (dopTdamp && dip.system == 0 && dip.MEtype == 0) 
          wt *= pT2damp / (dip.pT2 + pT2damp);
      }
    }

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < rndmPtr->flat());

}

//--------------------------------------------------------------------------

// Evolve a QED dipole end, either charged or photon. 

void TimeShower::pT2nextQED(double pT2begDip, double pT2sel, 
  TimeDipoleEnd& dip, Event& event) { 

  // Lower cut for evolution. Return if no evolution range.
  double pT2chgCut = (dip.chgType != 0 && abs(dip.chgType) != 3)
    ? pT2chgQCut : pT2chgLCut;
  double pT2endDip = max( pT2sel, pT2chgCut ); 
  if (pT2begDip < pT2endDip) return;   

  // Emission of photon or photon branching.
  bool hasCharge = (dip.chgType != 0);

  // Default values.
  double wtPSgam     = 0.;
  double chg2Sum     = 0.;
  double chg2SumL    = 0.;
  double chg2SumQ    = 0.;
  double zMinAbs     = 0.; 
  double emitCoefTot = 0.;

  // alpha_em at maximum scale provides upper estimate.
  double alphaEMmax  = alphaEM.alphaEM(renormMultFac * pT2begDip);
  double alphaEM2pi  = alphaEMmax / (2. * M_PI);

  // Emission: upper estimate for matrix element weighting; charge factor.
  if (hasCharge) {
    wtPSgam     = 2.;
    double chg2 = pow2(dip.chgType / 3.);

    // Determine overestimated z range. Find evolution coefficient.
    zMinAbs = 0.5 - sqrtpos( 0.25 - pT2endDip / dip.m2DipCorr );
    if (zMinAbs < SIMPLIFYROOT) zMinAbs = pT2endDip / dip.m2DipCorr;
    emitCoefTot = alphaEM2pi * chg2 * wtPSgam * log(1. / zMinAbs - 1.);

  // Branching: sum of squared charge factors for lepton and quark daughters.
  } else { 
    chg2SumL = max(0, min(3, nGammaToLepton));
    if      (nGammaToQuark > 4) chg2SumQ = 11. / 9.;
    else if (nGammaToQuark > 3) chg2SumQ = 10. / 9.;
    else if (nGammaToQuark > 2) chg2SumQ =  6. / 9.;
    else if (nGammaToQuark > 1) chg2SumQ =  5. / 9.;
    else if (nGammaToQuark > 0) chg2SumQ =  1. / 9.;

    // Total sum of squared charge factors. Find evolution coefficient. 
    chg2Sum     = chg2SumL + 3. * chg2SumQ; 
    emitCoefTot = alphaEM2pi * chg2Sum;
  }
  
  // Variables used inside evolution loop.
  dip.pT2 = pT2begDip;
  double wt; 
  
  // Begin evolution loop towards smaller pT values.
  do { 
 
    // Pick pT2 (in overestimated z range).
    dip.pT2 = dip.pT2 * pow(rndmPtr->flat(), 1. / emitCoefTot);
    wt = 0.;

    // Abort evolution if below cutoff scale, or below another branching.
    if ( dip.pT2 < pT2endDip) { dip.pT2 = 0.; return; }

    // Pick z according to dz/(1-z) or flat.
    if (hasCharge) dip.z = 1. - zMinAbs 
      * pow( 1. / zMinAbs - 1., rndmPtr->flat() );
    else           dip.z = rndmPtr->flat();
  
    // Do not accept branching if outside allowed z range.
    double zMin = 0.5 - sqrtpos( 0.25 - dip.pT2 / dip.m2DipCorr ); 
    if (zMin < SIMPLIFYROOT) zMin = dip.pT2 / dip.m2DipCorr;
    dip.m2 = dip.m2Rad + dip.pT2 / (dip.z * (1. - dip.z));
    if (dip.z > zMin && dip.z < 1. - zMin 
      && dip.m2 * dip.m2Dip < dip.z * (1. - dip.z) 
        * pow2(dip.m2Dip + dip.m2 - dip.m2Rec) 
      // For gamma -> f fbar also impose maximum mass.
      && (hasCharge || dip.m2 < m2MaxGamma) ) {

      // Photon emission: unique flavour choice.
      if (hasCharge) {
        dip.flavour = 22;
        dip.mFlavour = 0.;

      // Photon branching: either lepton or quark flavour choice.   
      } else {
        if (rndmPtr->flat() * chg2Sum < chg2SumL)  
          dip.flavour  = 9 + 2 * min(3, 1 + int(chg2SumL * rndmPtr->flat()));
        else { 
          double rndmQ = 9. * chg2SumQ * rndmPtr->flat();
          if      (rndmQ <  1.) dip.flavour = 1;
          else if (rndmQ <  5.) dip.flavour = 2;
          else if (rndmQ <  6.) dip.flavour = 3;
          else if (rndmQ < 10.) dip.flavour = 4;
          else                  dip.flavour = 5;
        }
        dip.mFlavour = particleDataPtr->m0(dip.flavour);
      }                      

      // No z weight, except threshold, if to do ME corrections later on.
      if (dip.MEtype > 0) { 
        wt = 1.;
        if (dip.flavour < 20 && dip.m2 < THRESHM2 * pow2(dip.mFlavour))
          wt = 0.; 

      // z weight for X -> X gamma.
      } else if (hasCharge) {
        wt = (1. + pow2(dip.z)) / wtPSgam;

      // z weight for gamma -> f fbar.
      } else {
        double beta  = sqrtpos( 1. - 4. * pow2(dip.mFlavour) / dip.m2 );
        wt = beta * ( pow2(dip.z) + pow2(1. - dip.z) );
      }

      // Correct to current value of alpha_EM.
      double alphaEMnow = alphaEM.alphaEM(renormMultFac * dip.pT2);
      wt *= (alphaEMnow / alphaEMmax);

      // Suppression factors for dipole to beam remnant.
      if (dip.isrType != 0 && useLocalRecoilNow) {
        BeamParticle& beam = (dip.isrType == 1) ? *beamAPtr : *beamBPtr;
        int iSys    = dip.system;
        double xOld = beam[iSys].x();
        double xNew = xOld * (1. + (dip.m2 - dip.m2Rad) / 
          (dip.m2Dip - dip.m2Rad));
        double xMaxAbs = beam.xMax(iSys);
        if (xMaxAbs < 0.) {
          infoPtr->errorMsg("Warning in TimeShower::pT2nextQED: "
          "xMaxAbs negative"); 
          return;
        }
 
        // Firstly reduce by PDF ratio.
        if (xNew > xMaxAbs) wt = 0.;
        else {
          int idRec     = event[dip.iRecoiler].id();
          double pdfOld = max ( TINYPDF, 
            beam.xfISR( iSys, idRec, xOld, factorMultFac * dip.pT2) ); 
          double pdfNew = 
            beam.xfISR( iSys, idRec, xNew, factorMultFac * dip.pT2); 
          wt *= min( 1., pdfNew / pdfOld); 
        }

        // Secondly optionally reduce by 4 pT2_hard / (4 pT2_hard + m2).
        if (dampenBeamRecoil) {
          double pT24 = 4. * event[dip.iRadiator].pT2();
          wt *= pT24 / (pT24 + dip.m2);
        }
      }

      // Optional dampening of large pT values in hard system.
      if (dopTdamp && dip.system == 0 && dip.MEtype == 0) 
        wt *= pT2damp / (dip.pT2 + pT2damp);
    }

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < rndmPtr->flat());  

}

//--------------------------------------------------------------------------

// Evolve a Hidden Valley dipole end. 

void TimeShower::pT2nextHV(double pT2begDip, double pT2sel, 
  TimeDipoleEnd& dip, Event& ) { 

  // Lower cut for evolution. Return if no evolution range.
  double pT2endDip = max( pT2sel, pT2hvCut ); 
  if (pT2begDip < pT2endDip) return;   

  // C_F * alpha_HV/2 pi.
  int    colvTypeAbs = abs(dip.colvType);
  double colvFac     = (colvTypeAbs == 1) ? CFHV : 0.5 * nCHV;
  double alphaHV2pi  = colvFac * (alphaHVfix / (2. * M_PI));

  // Determine overestimated z range. Find evolution coefficient.
  double zMinAbs = 0.5 - sqrtpos( 0.25 - pT2endDip / dip.m2DipCorr );
  if (zMinAbs < SIMPLIFYROOT) zMinAbs = pT2endDip / dip.m2DipCorr;
  double emitCoefTot = alphaHV2pi * 2. * log(1. / zMinAbs - 1.);
  
  // Variables used inside evolution loop.
  dip.pT2 = pT2begDip;
  double wt; 
  
  // Begin evolution loop towards smaller pT values.
  do { 
 
    // Pick pT2 (in overestimated z range).
    dip.pT2 = dip.pT2 * pow(rndmPtr->flat(), 1. / emitCoefTot);
    wt = 0.;

    // Abort evolution if below cutoff scale, or below another branching.
    if ( dip.pT2 < pT2endDip) { dip.pT2 = 0.; return; }

    // Pick z according to dz/(1-z).
    dip.z = 1. - zMinAbs * pow( 1. / zMinAbs - 1., rndmPtr->flat() );
  
    // Do not accept branching if outside allowed z range.
    double zMin = 0.5 - sqrtpos( 0.25 - dip.pT2 / dip.m2DipCorr ); 
    if (zMin < SIMPLIFYROOT) zMin = dip.pT2 / dip.m2DipCorr;
    dip.m2 = dip.m2Rad + dip.pT2 / (dip.z * (1. - dip.z));
    if (dip.z > zMin && dip.z < 1. - zMin 
      && dip.m2 * dip.m2Dip < dip.z * (1. - dip.z) 
        * pow2(dip.m2Dip + dip.m2 - dip.m2Rec) ) {

      // HV gamma or gluon emission: unique flavour choice.
      dip.flavour  = idHV;
      dip.mFlavour = mHV;

      // No z weight, except threshold, if to do ME corrections later on.
      if (dip.MEtype > 0) wt = 1.;

      // z weight for X -> X g_HV.
      else if (colvTypeAbs == 1) wt = (1. + pow2(dip.z)) / 2.;
      else wt = (1. + pow3(dip.z)) / 2.;
    }

    // Optional dampening of large pT values in hard system.
    if (dopTdamp && dip.system == 0 && dip.MEtype == 0) 
      wt *= pT2damp / (dip.pT2 + pT2damp);

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < rndmPtr->flat());  

}

//--------------------------------------------------------------------------

// ME corrections and kinematics that may give failure.
// Notation: radBef, recBef = radiator, recoiler before emission,
//           rad, rec, emt = radiator, recoiler, emitted efter emission.
//           (rad, emt distinguished by colour flow for g -> q qbar.) 

bool TimeShower::branch( Event& event, bool isInterleaved) {

  // Check if global recoil should be used.
  useLocalRecoilNow = !(globalRecoil && dipSel->system == 0 
    && partonSystemsPtr->sizeOut(0) <= nMaxGlobalRecoil);

  // Check if the first emission shoild be checked for removal
  bool canMergeFirst = (mergingHooksPtr != 0)
                     ? mergingHooksPtr->canVetoEmission() : false;

  // Find initial radiator and recoiler particles in dipole branching.
  int iRadBef      = dipSel->iRadiator;
  int iRecBef      = dipSel->iRecoiler;
  Particle& radBef = event[iRadBef]; 
  Particle& recBef = event[iRecBef];

  // Find their momenta, with special sum for global recoil.
  Vec4 pRadBef     = event[iRadBef].p(); 
  Vec4 pRecBef; 
  vector<int> iGRecBef, iGRec; 
  if (useLocalRecoilNow) pRecBef =  event[iRecBef].p(); 
  else {
    for (int i = 0; i < partonSystemsPtr->sizeOut( dipSel->system); ++i) { 
      int iG = partonSystemsPtr->getOut( dipSel->system, i);
      if (iG !=  dipSel->iRadiator) {
        iGRecBef.push_back(iG);
        pRecBef += event[iG].p();
      }
    }
  }

  // Default flavours and colour tags for new particles in dipole branching. 
  int idRad        = radBef.id();
  int idEmt        = dipSel->flavour; 
  int colRad       = radBef.col();
  int acolRad      = radBef.acol();
  int colEmt       = 0;
  int acolEmt      = 0;
  iSysSel          = dipSel->system;
  int iSysSelRec   = dipSel->systemRec;

  // Default OK for photon, photon_HV or gluon_HV emission.
  if (dipSel->flavour == 22 || dipSel->flavour == idHV) { 
  // New colour tag required for gluon emission.
  } else if (dipSel->flavour == 21 && dipSel->colType > 0) { 
    colEmt  = colRad;  
    colRad  = event.nextColTag();   
    acolEmt = colRad;
  } else if (dipSel->flavour == 21) { 
    acolEmt = acolRad;  
    acolRad = event.nextColTag();   
    colEmt  = acolRad;
  // New flavours for g -> q qbar; split colours.
  } else if (dipSel->colType > 0) {
    idEmt   = dipSel->flavour ;
    idRad   = -idEmt;
    colEmt  = colRad;
    colRad  = 0; 
  } else if (dipSel->colType < 0) {
    idEmt   = -dipSel->flavour ;
    idRad   = -idEmt;
    acolEmt = acolRad;
    acolRad = 0; 
  // New flavours for gamma -> f fbar, and maybe also colours.
  } else if (dipSel->gamType == 1 && rndmPtr->flat() > 0.5) {   
    idEmt   = -dipSel->flavour ;
    idRad   = -idEmt;
    if (idRad < 10) colRad = event.nextColTag(); 
    acolEmt = colRad;
  } else if (dipSel->gamType == 1) {   
    idEmt   = dipSel->flavour ;
    idRad   = -idEmt;
    if (idEmt < 10) colEmt = event.nextColTag(); 
    acolRad = colEmt;
  }

  // Construct kinematics in dipole rest frame: 
  // begin simple (like g -> g g).
  double pTorig       = sqrt( dipSel->pT2);
  double eRadPlusEmt  = 0.5 * (dipSel->m2Dip + dipSel->m2 - dipSel->m2Rec) 
    / dipSel->mDip;
  double e2RadPlusEmt = pow2(eRadPlusEmt);
  double pzRadPlusEmt = 0.5 * sqrtpos( pow2(dipSel->m2Dip - dipSel->m2 
    - dipSel->m2Rec) - 4. * dipSel->m2 * dipSel->m2Rec ) / dipSel->mDip;
  double pT2corr = dipSel->m2 * (e2RadPlusEmt * dipSel->z * (1. - dipSel->z)
                      - 0.25 * dipSel->m2) / pow2(pzRadPlusEmt);
  double pTcorr       = sqrtpos( pT2corr );
  double pzRad        = (e2RadPlusEmt * dipSel->z - 0.5 * dipSel->m2) 
                      / pzRadPlusEmt;
  double pzEmt        = (e2RadPlusEmt * (1. - dipSel->z) - 0.5 * dipSel->m2) 
                      / pzRadPlusEmt;
  double mRad         = dipSel->mRad;
  double mEmt         = 0.;

  // Kinematics reduction for q -> q gamma_v when m_q > 0 and m_gamma_v > 0.
  if ( abs(dipSel->colvType) == 1 && dipSel->mFlavour > 0.) {  
    mEmt              = dipSel->mFlavour;
    if (pow2(mRad + mEmt) > dipSel->m2) return false;
    double m2Emt      = pow2(mEmt);
    double lambda     = sqrtpos( pow2(dipSel->m2 - dipSel->m2Rad - m2Emt)
                      - 4. * dipSel->m2Rad * m2Emt );
    kRad              = 0.5 * (dipSel->m2 - lambda + m2Emt - dipSel->m2Rad) 
                      / dipSel->m2;
    kEmt              = 0.5 * (dipSel->m2 - lambda + dipSel->m2Rad - m2Emt)
                      / dipSel->m2; 
    pTorig           *= 1. - kRad - kEmt;
    pTcorr           *= 1. - kRad - kEmt;
    double pzMove     = kRad * pzRad - kEmt * pzEmt;
    pzRad            -= pzMove;
    pzEmt            += pzMove; 

  // Kinematics reduction for q -> q g/gamma/g_HV when m_q > 0. 
  } else if (abs(dipSel->colType) == 1 || dipSel->chgType != 0 
    || abs(dipSel->colvType) == 1) { 
    pTorig           *= 1. - dipSel->m2Rad / dipSel->m2; 
    pTcorr           *= 1. - dipSel->m2Rad / dipSel->m2; 
    pzRad            += pzEmt * dipSel->m2Rad / dipSel->m2;
    pzEmt            *= 1. - dipSel->m2Rad / dipSel->m2; 
 
  // Kinematics reduction for g -> q qbar or gamma -> f fbar when m_f > 0;
  } else if (abs(dipSel->flavour) < 20) {
    mEmt              = dipSel->mFlavour;
    mRad              = mEmt;
    double beta       = sqrtpos( 1. - 4. * pow2(mEmt) / dipSel->m2 );   
    pTorig           *= beta;
    pTcorr           *= beta;
    pzRad             = 0.5 * ( (1. + beta) * pzRad + (1. - beta) * pzEmt );
    pzEmt             = pzRadPlusEmt - pzRad;
  } 

  // Reject g emission where mass effects have reduced pT below cutoff.
  if (idEmt == 21 && pTorig < pTcolCut) return false;

  // Find rest frame and angles of original dipole.
  RotBstMatrix M;
  M.fromCMframe(pRadBef, pRecBef);

  // Evaluate coefficient of azimuthal asymmetry from gluon polarization.
  findAsymPol( event, dipSel);

  // Begin construction of new dipole kinematics: pick azimuthal angle.
  Vec4 pRad, pEmt, pRec;
  double wtPhi = 1.;
  do { 
    double phi = 2. * M_PI * rndmPtr->flat();

    // Define kinematics of branching in dipole rest frame.
    pRad = Vec4( pTcorr * cos(phi), pTcorr * sin(phi), pzRad, 
      sqrt( pow2(pTcorr) + pow2(pzRad) + pow2(mRad) ) );
    pEmt = Vec4( -pRad.px(), -pRad.py(), pzEmt,
      sqrt( pow2(pTcorr) + pow2(pzEmt) + pow2(mEmt) ) );
    pRec = Vec4( 0., 0., -pzRadPlusEmt, sqrt( pow2(pzRadPlusEmt) 
      + dipSel->m2Rec ) );

    // Rotate and boost dipole products to the event frame.
    pRad.rotbst(M);
    pEmt.rotbst(M);
    pRec.rotbst(M);

    // Azimuthal phi weighting: loop to new phi value if required.
    if (dipSel->asymPol != 0.) {
      Vec4 pAunt = event[dipSel->iAunt].p();
      double cosPhi = cosphi( pRad, pAunt, pRadBef );
      wtPhi = ( 1. + dipSel->asymPol * (2. * pow2(cosPhi) - 1.) )
        / ( 1. + abs(dipSel->asymPol) );
    } 
  } while (wtPhi < rndmPtr->flat()) ;

  // Kinematics when recoiler is initial-state parton.
  int isrTypeNow  = dipSel->isrType;
  int isrTypeSave = isrTypeNow;
  if (!useLocalRecoilNow) isrTypeNow = 0;
  if (isrTypeNow != 0) pRec = 2. * recBef.p() - pRec;

  // PS dec 2010: check if radiator has flexible normalization 
  bool isFlexible = dipSel->isFlexible;

  // Define new particles from dipole branching.
  double pTsel = sqrt(dipSel->pT2);
  Particle rad = Particle(idRad, 51, iRadBef, 0, 0, 0, 
    colRad, acolRad, pRad, mRad, pTsel); 
  Particle emt = Particle(idEmt, 51, iRadBef, 0, 0, 0,
    colEmt, acolEmt, pEmt, mEmt, pTsel);

  // Recoiler either in final or in initial state
  Particle rec = (isrTypeNow == 0)
    ? Particle(recBef.id(),  52, iRecBef, iRecBef, 0, 0, 
      recBef.col(), recBef.acol(), pRec, dipSel->mRec, pTsel) 
    : Particle(recBef.id(), -53, 0, 0, iRecBef, iRecBef, 
      recBef.col(), recBef.acol(), pRec, 0., 0.); 

  // ME corrections can lead to branching being rejected.
  if (dipSel->MEtype > 0) {
    Particle& partner = (dipSel->iMEpartner == iRecBef) 
      ? rec : event[dipSel->iMEpartner];
    if ( findMEcorr( dipSel, rad, partner, emt) < rndmPtr->flat() )
      return false;
  }

  // Rescatter: if the recoiling partner is not in the same system
  //            as the radiator, fix up intermediate systems (can lead
  //            to emissions being vetoed)
  if (allowRescatter && FIXRESCATTER && isInterleaved 
    && iSysSel != iSysSelRec) {
    Vec4 pNew = rad.p() + emt.p();
    if (!rescatterPropagateRecoil(event, pNew)) return false;
  }

  // Save properties to be restored in case of user-hook veto of emission.
  int eventSizeOld = event.size();
  int iRadStatusV  = event[iRadBef].status();
  int iRadDau1V    = event[iRadBef].daughter1();
  int iRadDau2V    = event[iRadBef].daughter2();
  int iRecStatusV  = event[iRecBef].status();
  int iRecMot1V    = event[iRecBef].mother1();
  int iRecMot2V    = event[iRecBef].mother2();
  int iRecDau1V    = event[iRecBef].daughter1();
  int iRecDau2V    = event[iRecBef].daughter2();
  int beamOff1     = 1 + beamOffset;
  int beamOff2     = 2 + beamOffset;
  int ev1Dau1V     = event[beamOff1].daughter1();
  int ev2Dau1V     = event[beamOff2].daughter1();

  // Shower may occur at a displaced vertex.
  if (radBef.hasVertex()) {
    rad.vProd( radBef.vProd() );
    emt.vProd( radBef.vProd() );
  }
  if (recBef.hasVertex()) rec.vProd( recBef.vProd() );

  // Put new particles into the event record.
  // Mark original dipole partons as branched and set daughters/mothers.
  int iRad = event.append(rad);
  int iEmt = event.append(emt);
  event[iRadBef].statusNeg();
  event[iRadBef].daughters( iRad, iEmt); 
  int iRec = 0; 
  if (useLocalRecoilNow) {
    iRec = event.append(rec);
    if (isrTypeNow == 0) {
      event[iRecBef].statusNeg();
      event[iRecBef].daughters( iRec, iRec);
    } else {
      event[iRecBef].mothers( iRec, iRec);
      event[iRec].mothers( iRecMot1V, iRecMot2V);  
      if (iRecMot1V == beamOff1) event[beamOff1].daughter1( iRec);  
      if (iRecMot1V == beamOff2) event[beamOff2].daughter1( iRec);
    } 
 
  // Global recoil: need to find relevant rotation+boost for recoilers:
  // boost+rotate to rest frame, boost along z axis, rotate+boost back.
  } else {
    RotBstMatrix MG = M;
    MG.invert();
    double pzRecBef = -0.5 * sqrtpos( pow2(dipSel->m2Dip - dipSel->m2Rad 
      - dipSel->m2Rec) - 4. * dipSel->m2Rad * dipSel->m2Rec ) / dipSel->mDip;
    double eRecBef  = sqrt( pow2(pzRecBef) + dipSel->m2Rec);
    double pzRecAft = -pzRadPlusEmt; 
    double eRecAft  = sqrt( pow2(pzRecAft) + dipSel->m2Rec);
    MG.bst( Vec4(0., 0., pzRecBef, eRecBef), Vec4(0., 0., pzRecAft, eRecAft) );
    MG.rotbst( M);

    // Global recoil: copy particles, and rotate+boost momenta (not vertices).
    for (int iG = 0; iG < int(iGRecBef.size()); ++iG) {
      iRec = event.copy( iGRecBef[iG], 52);
      iGRec.push_back( iRec);
      Vec4 pGRec = event[iRec].p();
      pGRec.rotbst( MG);
      event[iRec].p( pGRec);
    } 
  }  

  // Allow veto of branching. If so restore event record to before emission.
  bool inResonance = (partonSystemsPtr->getInA(iSysSel) == 0) ? true : false;
  if ( (canVetoEmission && userHooksPtr->doVetoFSREmission( eventSizeOld, 
    event, iSysSel, inResonance))
    || (canMergeFirst && mergingHooksPtr->doVetoEmission( event )) ) {
    event.popBack( event.size() - eventSizeOld);
    event[iRadBef].status( iRadStatusV);
    event[iRadBef].daughters( iRadDau1V, iRadDau2V);
    if (useLocalRecoilNow && isrTypeNow == 0) {
      event[iRecBef].status( iRecStatusV);
      event[iRecBef].daughters( iRecDau1V, iRecDau2V);
    } else if (useLocalRecoilNow) {
      event[iRecBef].mothers( iRecMot1V, iRecMot2V);
      if (iRecMot1V == beamOff1) event[beamOff1].daughter1( ev1Dau1V);
      if (iRecMot1V == beamOff2) event[beamOff2].daughter1( ev2Dau1V);
    } else {
      for (int iG = 0; iG < int(iGRecBef.size()); ++iG) {
        event[iGRecBef[iG]].statusPos();
        event[iGRecBef[iG]].daughters( 0, 0);
      }
    }
    return false;
  }
 
  // For global recoil restore the one nominal recoiler, for bookkeeping.
  if (!useLocalRecoilNow) {
    iRec = iRecBef;
    for (int iG = 0; iG < int(iGRecBef.size()); ++iG) 
    if (iGRecBef[iG] == iRecBef) iRec = iGRec[iG];
  }

  // For initial-state recoiler also update beam and sHat info.
  if (isrTypeNow != 0) {
    BeamParticle& beamRec = (isrTypeNow == 1) ? *beamAPtr : *beamBPtr;
    double xOld = beamRec[iSysSelRec].x();
    double xRec = 2. * pRec.e() / (beamAPtr->e() + beamBPtr->e());
    beamRec[iSysSelRec].iPos( iRec);
    beamRec[iSysSelRec].x( xRec); 
    partonSystemsPtr->setSHat( iSysSelRec,
    partonSystemsPtr->getSHat(iSysSelRec) * xRec / xOld);
  }

  // Photon emission: update to new dipole ends; add new photon "dipole".
  if (dipSel->flavour == 22) { 
    dipSel->iRadiator = iRad;
    dipSel->iRecoiler = iRec;
    // When recoiler was uncharged particle, in resonance decays,
    // assign recoil to emitted photons. 
    if (recoilToColoured && inResonance && event[iRec].chargeType() == 0) 
      dipSel->iRecoiler = iEmt;
    dipSel->pTmax = pTsel;
    if (doQEDshowerByGamma) dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, 
      pTsel, 0, 0, 1, 0, iSysSel, 0));
 
  // Gluon emission: update both dipole ends and add two new ones.
  } else if (dipSel->flavour == 21) { 
    dipSel->iRadiator  = iRad;
    dipSel->iRecoiler  = iEmt;
    dipSel->systemRec  = iSysSel;
    dipSel->isrType    = 0;
    dipSel->pTmax      = pTsel;
    // Optionally also kill ME corrections after first emission.
    if (!doMEafterFirst) dipSel->MEtype = 0;
    // PS dec 2010: check normalization of radiating dipole 
    // Dipole corresponding to the newly created color tag has normal strength
    double flexFactor  = (isFlexible) ? dipSel->flexFactor : 1.0;
    dipSel->isFlexible = false;
    for (int i = 0; i < int(dipEnd.size()); ++i) {
      if (dipEnd[i].iRadiator == iRecBef && dipEnd[i].iRecoiler == iRadBef 
        && dipEnd[i].colType != 0) {
        dipEnd[i].iRadiator = iRec;
        dipEnd[i].iRecoiler = iEmt;
        // Optionally also kill ME corrections after first emission.
        if (!doMEafterFirst) dipEnd[i].MEtype = 0;
        // Strive to match colour to anticolour inside closed system.
        if ( !isFlexible && dipEnd[i].colType * dipSel->colType > 0) 
          dipEnd[i].iRecoiler = iRad;
        dipEnd[i].pTmax = pTsel;	
	// PS dec 2010: if the (iRadBef,iRecBef) dipole was flexible, the
	// same should be true for this (opposite) end. If so, this end keeps 
	// the modified normalization, so we shouldn't need to do anything. 
      }
    }
    int colType = (dipSel->colType > 0) ? 2 : -2 ;
    // When recoiler was uncoloured particle, in resonance decays, 
    // assign recoil to coloured particle.
    int iRecMod = iRec;
    if (recoilToColoured && inResonance && event[iRec].col() == 0 
      && event[iRec].acol() == 0) iRecMod = iRad;
    dipEnd.push_back( TimeDipoleEnd(iEmt, iRecMod, pTsel,  
       colType, 0, 0, isrTypeSave, iSysSel, 0));
    dipEnd.back().systemRec = iSysSelRec;
    // PS dec 2010: the (iEmt,iRec) dipole "inherits" flexible normalization
    if (isFlexible) {
      dipEnd.back().isFlexible = true;
      dipEnd.back().flexFactor = flexFactor;
    }
    dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel, 
      -colType, 0, 0, 0, iSysSel, 0));
    
  // Gluon branching to q qbar: update current dipole and other of gluon.
  } else if (dipSel->colType != 0) {
    for (int i = 0; i < int(dipEnd.size()); ++i) {
      // Strive to match colour to anticolour inside closed system.
      if ( !isFlexible && dipEnd[i].iRecoiler == iRadBef 
        && dipEnd[i].colType * dipSel->colType < 0 ) 
        dipEnd[i].iRecoiler = iEmt;
      if (dipEnd[i].iRadiator == iRadBef && abs(dipEnd[i].colType) == 2) {
        dipEnd[i].colType /= 2;
        if (dipEnd[i].system != dipEnd[i].systemRec) continue;

        // Note: gluino -> quark + squark gives a deeper radiation dip than
        // the more obvious alternative photon decay, so is more realistic.
        dipEnd[i].MEtype = 66;
        if (&dipEnd[i] == dipSel) dipEnd[i].iMEpartner = iRad;
        else                      dipEnd[i].iMEpartner = iEmt;
      }
    }
    dipSel->iRadiator = iEmt;
    dipSel->iRecoiler = iRec;
    dipSel->pTmax     = pTsel;
    
    // Gluon branching to q qbar: also add two charge dipole ends.
    // Note: gluino -> quark + squark gives a deeper radiation dip than
    // the more obvious alternative photon decay, so is more realistic.
    if (doQEDshowerByQ) {
      int chgType = event[iRad].chargeType(); 
      dipEnd.push_back( TimeDipoleEnd(iRad, iEmt, pTsel, 
        0,  chgType, 0, 0, iSysSel, 66, iEmt));
      dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel, 
        0, -chgType, 0, 0, iSysSel, 66, iRad));
    }

  // Photon branching to f fbar: inactivate photon "dipole";
  // optionally add new charge and colour dipole ends. 
  } else if (dipSel->gamType != 0) {
    dipSel->gamType = 0;
    int chgType = event[iRad].chargeType(); 
    int colType = event[iRad].colType();
    // MEtype = 102 for charge in vector decay.
    if ( chgType != 0 && ( ( doQEDshowerByQ && colType != 0 )  
      || ( doQEDshowerByL && colType == 0 ) ) ) { 
      dipEnd.push_back( TimeDipoleEnd(iRad, iEmt, pTsel, 
        0,  chgType, 0, 0, iSysSel, 102, iEmt));
      dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel, 
        0, -chgType, 0, 0, iSysSel, 102, iRad));
    }
    // MEtype = 11 for colour in vector decay.
    if (colType != 0 && doQCDshower) {
      dipEnd.push_back( TimeDipoleEnd(iRad, iEmt, pTsel, 
         colType, 0, 0, 0, iSysSel, 11, iEmt));
      dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel, 
        -colType, 0, 0, 0, iSysSel, 11, iRad));
    }

  // Photon_HV emission: update to new dipole ends.
  } else if (dipSel->flavour == 4900022) { 
    dipSel->iRadiator = iRad;
    dipSel->iRecoiler = iRec;
    dipSel->pTmax = pTsel;

  // Gluon_HV emission: update to new dipole ends.
  } else if (dipSel->flavour == 4900021) { 
    dipSel->iRadiator = iRad;
    dipSel->iRecoiler = iEmt;
    dipSel->pTmax     = pTsel;
    for (int i = 0; i < int(dipEnd.size()); ++i) 
    if (dipEnd[i].iRadiator == iRecBef && dipEnd[i].iRecoiler == iRadBef 
      && dipEnd[i].isHiddenValley) {
      dipEnd[i].iRadiator = iRec;
      dipEnd[i].iRecoiler = iEmt;
      dipEnd[i].pTmax     = pTsel;
    }
    int colvType = (dipSel->colvType > 0) ? 2 : -2 ;
    dipEnd.push_back( TimeDipoleEnd(iEmt, iRec, pTsel,  
       0, 0, 0, isrTypeSave, iSysSel, 0, -1, false, true, colvType) );
    dipEnd.back().systemRec = iSysSelRec;
    dipEnd.push_back( TimeDipoleEnd(iEmt, iRad, pTsel, 
      0, 0, 0, 0, iSysSel, 0, -1, false, true, -colvType) );
  }

  // Copy or set lifetime for new final state. 
  if (event[iRad].id() == event[iRadBef].id()) {
    event[iRad].tau( event[iRadBef].tau() );
  } else {
    event[iRad].tau( event[iRad].tau0() * rndmPtr->exp() );
    event[iEmt].tau( event[iEmt].tau0() * rndmPtr->exp() );
  } 
  event[iRec].tau( event[iRecBef].tau() );

  // Now update other dipoles that also involved the radiator or recoiler.
  for (int i = 0; i < int(dipEnd.size()); ++i) {
    // PS dec 2010: if radiator was flexible and now is normal, there may
    // be other flexible dipoles that need updating.
    if (isFlexible && !dipSel->isFlexible && dipEnd[i].isFlexible) {
      if (dipEnd[i].iRecoiler  == iRadBef) dipEnd[i].iRecoiler = iEmt;
      if (dipEnd[i].iRadiator  == iRadBef) {
	dipEnd[i].iRadiator = iEmt;
	if (dipEnd[i].colType == 1 && dipSel->flavour == 21) 
          dipEnd[i].colType = 2;
	if (dipEnd[i].colType ==-1 && dipSel->flavour == 21) 
          dipEnd[i].colType =-2;
      }
    }
    if (dipEnd[i].iRadiator  == iRadBef) dipEnd[i].iRadiator  = iRad;
    if (dipEnd[i].iRecoiler  == iRadBef) dipEnd[i].iRecoiler  = iRad;
    if (dipEnd[i].iMEpartner == iRadBef) dipEnd[i].iMEpartner = iRad;
    if (useLocalRecoilNow) {
      if (dipEnd[i].iRadiator  == iRecBef) dipEnd[i].iRadiator  = iRec;
      if (dipEnd[i].iRecoiler  == iRecBef) dipEnd[i].iRecoiler  = iRec;
      if (dipEnd[i].iMEpartner == iRecBef) dipEnd[i].iMEpartner = iRec;
    } else {
      for (int iG = 0; iG < int(iGRecBef.size()); ++iG) {
        if (dipEnd[i].iRadiator  == iGRecBef[iG]) 
            dipEnd[i].iRadiator  =  iGRec[iG];
        if (dipEnd[i].iRecoiler  == iGRecBef[iG]) 
            dipEnd[i].iRecoiler  =  iGRec[iG];
        if (dipEnd[i].iMEpartner == iGRecBef[iG]) 
            dipEnd[i].iMEpartner =  iGRec[iG];
      }
    }
  }

  // PS Apr 2011
  // Update any junctions downstream of this branching (if necessary)
  // (This happens, e.g., via LHEF, when adding showers to intermediate 
  //  coloured resonances whose decays involved junctions)
  for (int iJun = 0; iJun < event.sizeJunction(); iJun++) {
    // Number of incoming colour lines for this junction.
    int nIncoming = (event.kindJunction(iJun)-1)/2;
    // Check radiator colour or anticolour, depending on junction kind
    // (if junction, incoming = anticolours, and vice versa)
    int colChk = 0; 
    colChk = ( event.kindJunction(iJun) % 2 == 0 )
           ? event[iRadBef].col() : event[iRadBef].acol();
    // Loop over incoming junction ends
    for (int iCol = 0; iCol < nIncoming; iCol++) {      
      int colJun = event.colJunction( iJun, iCol);      
      // If match, update junction end with new upstream (anti)colour 
      if (colJun == colChk) {
	int colNew = 0;
	if ( event.kindJunction(iJun) % 2 == 0 ) colNew = colRad;
	else colNew = acolRad;
	event.colJunction( iJun, iCol, colNew );
      }
    }
  }

  // Finally update the list of all partons in all systems.
  partonSystemsPtr->replace(iSysSel, iRadBef, iRad);  
  partonSystemsPtr->addOut(iSysSel, iEmt);
  if (useLocalRecoilNow) 
    partonSystemsPtr->replace(iSysSelRec, iRecBef, iRec);
  else {
    for (int iG = 0; iG < int(iGRecBef.size()); ++iG) 
    partonSystemsPtr->replace(iSysSel, iGRecBef[iG], iGRec[iG]);
  }

  // Done. 
  return true;

}

//--------------------------------------------------------------------------

// Rescatter: If a dipole stretches between two different systems, those
//            systems will no longer locally conserve momentum. These
//            imbalances become problematic when ISR or primordial kT
//            is switched on as these steps involve Lorentz boosts.
//
//            'rescatterPropagateRecoil' tries to fix momentum in all
//            systems by propogating recoil momentum through all
//            intermediate systems. As the momentum transfer is already
//            defined, this can lead to internal lines gaining a
//            virtuality.

// Useful definitions for a pair of integers and a vector of pairs
typedef pair < int, int >  pairInt;
typedef vector < pairInt > vectorPairInt;

//--------------------------------------------------------------------------

// findParentSystems
//  Utility routine to find all parent systems of a given system
//  Returns a vector of pairs of integers with:
//   a) The system index, including the starting system (negative
//      if (b) points to a parent system, positive if (b) points
//      to a daughter system
//   b) The event record index that is the path out of the system
//      (if forwards == false, this is an incoming parton to the
//      system, and is +ve if side A or -ve if side B,
//      if forwards == true, this is an outgoing parton from the
//      system).
//  Returns as empty vector on failure
//  Note: this assumes single rescattering only and therefore only
//        one possible parent system

inline vectorPairInt findParentSystems(const int sys, 
  Event& event, PartonSystems* partonSystemsPtr, bool forwards) {

  vectorPairInt parentSystems;
  parentSystems.reserve(10);

  int iSysCur = sys;
  while (true) {
    // Get two incoming partons
    int iInA = partonSystemsPtr->getInA(iSysCur);
    int iInB = partonSystemsPtr->getInB(iSysCur);

    // Check if either of these links to another system
    int iIn = 0;
    if (event[iInA].isRescatteredIncoming()) iIn =  iInA;
    if (event[iInB].isRescatteredIncoming()) iIn = -iInB;

    // Save the current system to the vector
    parentSystems.push_back( pairInt(-iSysCur, iIn) ); 
    if (iIn == 0) break;

    int iInAbs  = abs(iIn);
    int iMother = event[iInAbs].mother1();
    iSysCur     = partonSystemsPtr->getSystemOf(iMother);
    if (iSysCur == -1) {
      parentSystems.clear();
      break;
    }
  } // while (true)

  // If forwards is set, change all event record indices to go to daughter
  // systems rather than parent systems
  if (forwards) {
    vectorPairInt::reverse_iterator rit;
    for (rit = parentSystems.rbegin(); rit < (parentSystems.rend() - 1);
         ++rit) {
      pairInt &cur  = *rit;
      pairInt &next = *(rit + 1);
      cur.first     = -cur.first;
      cur.second    = (next.second < 0) ? -event[abs(next.second)].mother1() :
                                           event[abs(next.second)].mother1();
    }
  } 

  return parentSystems;
}

//--------------------------------------------------------------------------

// rescatterPropagateRecoil
//  Fix up momentum in all intermediate systems when radiator and recoiler
//  systems are different. The strategy is to look at all parent systems
//  from the radiator system and the recoiler system and find where they
//  intersect.

bool TimeShower::rescatterPropagateRecoil( Event& event, Vec4& pNew) {

  // Some useful variables for later
  int  iRadBef    = dipSel->iRadiator;
       iSysSel    = dipSel->system;
  int  iSysSelRec = dipSel->systemRec;
  Vec4 pImbal     = pNew - event[iRadBef].p();

  // Store changes locally at first in case we veto the branching
  // eventMod stores index into the event record and new 4-vector
  vector < pair < int, Vec4 > > eventMod;
  eventMod.reserve(10);
  // systemMod stores system index (iSys) and system-parton index (iMem)
  //   iMem >=  0 - index into outgoing partons (iOut)
  //   iMem == -1 - incoming A
  //   iMem == -2 - incoming B
  vectorPairInt systemMod;
  systemMod.reserve(10);

  // Find all parent systems from radiating and recoiling systems
  vectorPairInt radParent = findParentSystems(iSysSel, event,
                                              partonSystemsPtr, false);
  vectorPairInt recParent = findParentSystems(iSysSelRec, event,
                                              partonSystemsPtr, true);
  if (radParent.size() == 0 || recParent.size() == 0) {
    // This should never happen
    infoPtr->errorMsg("Error in TimeShower::rescatterPropagateRecoil: "
      "couldn't find parent system; branching vetoed");
    return false;
  }
  // Find the system that connects radiating and recoiling system
  bool foundPath = false;
  unsigned int iRadP = 0; 
  unsigned int iRecP = 0;
  for (iRadP = 0; iRadP < radParent.size(); iRadP++) {
    for (iRecP = 0; iRecP < recParent.size(); iRecP++)
      if (abs(radParent[iRadP].first) == abs(recParent[iRecP].first)) {
        foundPath = true;
        break;
      }
    if (foundPath) break;
  }
  if (!foundPath) {
    // Can fail e.g. for QED dipoles where there is no connection
    // between radiator and recoiler systems
    infoPtr->errorMsg("Warning in TimeShower::rescatterPropagateRecoil: "
      "couldn't find recoil path; branching vetoed");
    return false;
  }

  // Join together to form complete path from radiating system
  // to recoiling system
  vectorPairInt path;
  if (radParent.size() > 1)
    path.assign(radParent.begin(), radParent.begin() + iRadP);
  if (recParent.size() > 1)
    path.insert(path.end(), recParent.rend() - iRecP - 1,
                recParent.rend() - 1);

  // Follow the path fixing up momenta as we go
  for (unsigned int i = 0; i < path.size(); i++) {
    // Line out of the current system
    bool isIncoming  = (path[i].first < 0) ? true : false;
    int  iSysCur     = abs(path[i].first);
    bool isIncomingA = (path[i].second > 0) ? true : false;
    int  iLink       = abs(path[i].second);

    int iMemCur;
    if (isIncoming) iMemCur = (isIncomingA) ? -1 : -2;
    else {
      iMemCur = -1;
      for (int j = 0; j < partonSystemsPtr->sizeOut(iSysCur); j++)
        if (partonSystemsPtr->getOut(iSysCur, j) == iLink) {
          iMemCur = j;
          break;
        }
      if (iMemCur == -1) {
        // This should never happen
        infoPtr->errorMsg("Error in TimeShower::rescatterPropagateRecoil: "
          "couldn't find parton system; branching vetoed");
        return false;
      }
    }

    Vec4 pMod = (isIncoming) ? event[iLink].p() + pImbal :
                               event[iLink].p() - pImbal;
    eventMod.push_back(pair <int, Vec4> (iLink, pMod));
    systemMod.push_back(pairInt(iSysCur, iMemCur));

    // Calculate sHat of iSysCur
    int  iInCurA = partonSystemsPtr->getInA(iSysCur);
    int  iInCurB = partonSystemsPtr->getInB(iSysCur);
    Vec4 pTotCur = event[iInCurA].p() + event[iInCurB].p();

    // If iMemCur is -1 or -2, then we must have changed the sHat of iSysCur
    if (iMemCur < 0) pTotCur += (isIncoming) ? pImbal : -pImbal;
    double sHatCur = pTotCur.m2Calc();
 
    // The fixed-up incoming and outgoing partons should not have
    // too large a virtuality in relation to the system mass-square. 
    if (abs(pMod.m2Calc()) > MAXVIRTUALITYFRACTION * sHatCur) {
      infoPtr->errorMsg("Warning in TimeShower::rescatterPropagateRecoil: "
        "virtuality much larger than sHat; branching vetoed");
      return false;
    }

    // Outgoing ones should also not have too large negative energy  
    // in the rest frame of the system.
    if (!isIncoming && pMod * pTotCur < -MAXNEGENERGYFRACTION * sHatCur) {
      infoPtr->errorMsg("Warning in TimeShower::rescatterPropagateRecoil: "
        "rest frame energy too negative; branching vetoed");
      return false;
    }

    // Veto negative sHat.
    if (sHatCur < 0.0) {
      infoPtr->errorMsg("Warning in TimeShower::rescatterPropagateRecoil: "
        "sHat became negative; branching vetoed");
      return false;
    }

    // Line into the new current system
    iLink   = (isIncoming) ? event[iLink].mother1()  :
                             event[iLink].daughter1();
    iSysCur = partonSystemsPtr->getSystemOf(iLink, true);

    if (!isIncoming) iMemCur = (isIncomingA) ? -1 : -2;
    else {
      iMemCur = -1;
      for (int j = 0; j < partonSystemsPtr->sizeOut(iSysCur); j++)
        if (partonSystemsPtr->getOut(iSysCur, j) == iLink) {
          iMemCur = j;
          break;
        }
      if (iMemCur == -1) {
        // This should never happen
        infoPtr->errorMsg("Error in TimeShower::rescatterPropagateRecoil: "
          "couldn't find parton system; branching vetoed");
        return false;
      }
    }

    pMod = (isIncoming) ? event[iLink].p() + pImbal :
                          event[iLink].p() - pImbal;
    eventMod.push_back(pair <int, Vec4> (iLink, pMod));
    systemMod.push_back(pairInt(iSysCur, iMemCur));

    // Calculate sHat of iSysCur
    iInCurA = partonSystemsPtr->getInA(iSysCur);
    iInCurB = partonSystemsPtr->getInB(iSysCur);
    pTotCur = event[iInCurA].p() + event[iInCurB].p();

    // If iMemCur is -1 or -2, then we must have changed the sHat of iSysCur
    if (iMemCur < 0) pTotCur += (isIncoming) ? pImbal : -pImbal;
    sHatCur = pTotCur.m2Calc();
 
    // The fixed-up incoming and outgoing partons should not have
    // too large a virtuality in relation to the system mass-square. 
    if (abs(pMod.m2Calc()) > MAXVIRTUALITYFRACTION * sHatCur) {
      infoPtr->errorMsg("Warning in TimeShower::rescatterPropagateRecoil: "
        "virtuality much larger than sHat; branching vetoed");
      return false;
    }

    // Outgoing ones should also not have too large negative energy
    // in the rest frame of the system.
    if (!isIncoming && pMod * pTotCur < -MAXNEGENERGYFRACTION * sHatCur) {
      infoPtr->errorMsg("Warning in TimeShower::rescatterPropagateRecoil: "
        "rest frame energy too negative; branching vetoed");
      return false;
    }

    // Veto negative sHat
    if (sHatCur < 0.0) {
      infoPtr->errorMsg("Warning in TimeShower::rescatterPropagateRecoil: "
        "sHat became negative; branching vetoed");
      return false;
    }

    // Do negative energy veto
    if (VETONEGENERGY && pMod.e() < 0.0) {
      infoPtr->errorMsg("Warning in TimeShower::rescatterPropagateRecoil: "
        "energy became negative; branching vetoed");
      return false;
    }
    
  } // for (unsigned int i = 0; i < path.size(); i++)

  // If no vetos by this point, apply the changes to the event record
  // An incoming particle with changed momentum is given status code -54,
  // an outgoing particle with changed momentum is given status code -55
  for (unsigned int i = 0; i < eventMod.size(); i++) {
    int idx    = eventMod[i].first;
    Vec4 &pMod = eventMod[i].second;
    int iSys   = systemMod[i].first;
    int iMem   = systemMod[i].second;

    // If incoming to a process then set the copy to be the mother
    if (event[idx].isRescatteredIncoming()) {
      int mother1 = event[idx].mother1();
      idx = event.copy(idx, -54);
      event[mother1].daughters(idx, idx);
    
      // Update beam information if necessary 
      double eCM = sqrt(m2( beamAPtr->p(), beamBPtr->p()));
      if        (iMem == -1) {
        partonSystemsPtr->setInA(iSys, idx);
        (*beamAPtr)[iSys].x((pMod.e() + pMod.pz()) / eCM);
        (*beamAPtr)[iSys].m(pMod.mCalc());
        (*beamAPtr)[iSys].p(pMod);
        (*beamAPtr)[iSys].iPos(idx);
      } else if (iMem == -2) {
        partonSystemsPtr->setInB(iSys, idx);
        (*beamBPtr)[iSys].x((pMod.e() - pMod.pz()) / eCM);
        (*beamBPtr)[iSys].m(pMod.mCalc());
        (*beamBPtr)[iSys].p(pMod);
        (*beamBPtr)[iSys].iPos(idx);
      } else {
        // This should never happen
        infoPtr->errorMsg("Error in TimeShower::rescatterPropagateRecoil: "
        "internal bookeeping error");
      }

    // Otherwise set the new event record entry to be the daughter
    } else {
      int daughter1 = event[idx].daughter1();
      idx = event.copy(idx, 55);
      event[idx].statusNeg();
      event[daughter1].mothers(idx, idx);

      partonSystemsPtr->setOut(iSys, iMem, idx);
    }

    event[idx].p( eventMod[i].second );
    event[idx].m( event[idx].mCalc() );
  }

  return true;
}


//--------------------------------------------------------------------------

// Find class of QCD ME correction.
// MEtype classification follow codes in Norrbin article,
// additionally -1 = try to find type, 0 = no ME corrections.
// Warning: not yet tried out to do a correct assignment in 
// arbitrary multiparton configurations! ??

void TimeShower::findMEtype( Event& event, TimeDipoleEnd& dip) {

  // Initial value. Mark if no ME corrections to be applied.
  bool setME = true;
  if (!doMEcorrections) setME = false; 
  int iMother  = event[dip.iRadiator].mother1();
  int iMother2 = event[dip.iRadiator].mother2();

  // Allow ME corrections for Hidden Valley pair in 2 -> 2.
  if (dip.isHiddenValley && event[dip.iRecoiler].id() 
    == -event[dip.iRadiator].id());

  // Else no ME corrections in 2 -> n processes.
  else {
    if (iMother2 != iMother && iMother2 != 0) setME = false;
    if (event[dip.iRecoiler].mother1() != iMother)  setME = false;    
    if (event[dip.iRecoiler].mother2() != iMother2) setME = false; 
  }   

  // No ME corrections for recoiler in initial state.
  if (event[dip.iRecoiler].status() < 0) setME = false;  

  // No ME corrections for recoiler not in same system
  if (dip.system != dip.systemRec) setME = false;

  // Done if no ME to be set.
  if (!setME) {
    dip.MEtype = 0;
    return;
  } 

  // If no ME partner set, assume it is the recoiler.
  if (dip.iMEpartner < 0) dip.iMEpartner = dip.iRecoiler;

  // Now begin processing of colour dipole, including Hidden Valley.
  if (dip.colType != 0 || dip.colvType != 0) {
    bool isHiddenColour = (dip.colvType != 0);

    // Find daughter types (may or may not be used later on).
    int idDau1      = event[dip.iRadiator].id();
    int idDau2      = event[dip.iMEpartner].id();
    int dau1Type    = findMEparticle(idDau1, isHiddenColour);
    int dau2Type    = findMEparticle(idDau2, isHiddenColour);
    int minDauType  = min(dau1Type, dau2Type);
    int maxDauType  = max(dau1Type, dau2Type);

    // Reorder dipole ends in kinematics. Split ME expression in two sides.
    dip.MEorder     = (dau2Type >= dau1Type);
    dip.MEsplit     = (maxDauType <= 6); 
    dip.MEgluinoRec = false;
 
    // If type already set (or set not to have) then done.
    if (minDauType == 0 && dip.MEtype < 0) dip.MEtype = 0;
    if (dip.MEtype >= 0) return;
    dip.MEtype = 0;

    // For H -> gg -> ggg we found that DGLAP kernels do better than eikonal.
    if (dau1Type == 4 && dau2Type == 4) return; 

    // Find mother type. 
    int idMother = 0;
    if ( event[dip.iRecoiler].mother1() == iMother && iMother >= 0) 
      idMother = event[iMother].id();
    int motherType = (idMother != 0) 
      ? findMEparticle(idMother, isHiddenColour) : 0;

    // When a mother if not known then use colour and spin content to guess.
    if (motherType == 0) {
      int col1  = event[dip.iRadiator].col();
      int acol1 = event[dip.iRadiator].acol();
      int col2  = event[dip.iMEpartner].col();
      int acol2 = event[dip.iMEpartner].acol();
      // spinT = 0/1 = integer or half-integer.
      int spinT = ( event[dip.iRadiator].spinType() 
                + event[dip.iMEpartner].spinType() )%2;
      // Colour singlet mother.
      if ( col1 == acol2 && acol1 == col2 ) 
        motherType = (spinT == 0) ? 7 : 9;
      // Colour octet mother.
      else if ( (col1 == acol2 && acol1 != 0 && col2 != 0)
        || (acol1 == col2 && col1 != 0 && acol2 != 0) )
        motherType = (spinT == 0) ? 4 : 5; 
      // Colour triplet mother.
      else if ( (col1 == acol2 && acol1 != col2)  
        || (acol1 == col2 && col1 != acol2) ) 
        motherType = (spinT == 0) ? 2 : 1;
      // If no colours are matched then cannot have common mother, so done.  
      else return;      
    }

    // Now start from default, which is eikonal ME corrections, 
    // and try to find matching ME cases below.
    int MEkind = 0;
    int MEcombi = 4;
    dip.MEmix = 0.5;

    // Hidden Valley with massive gamma_v covered by two special cases.
    if (isHiddenColour && brokenHVsym) {
      MEkind = (dau2Type == 0 || dau2Type > 6) ? 30 : 31;
      dip.MEtype = 5 * MEkind + 1; 
      return;
    }

    // Triplet recoiling against gluino needs enhanced radiation
    // to match to matrix elements.
    dip.MEgluinoRec = (dau1Type >= 1 && dau1Type <= 3 && dau2Type == 5);

    // Vector/axial vector -> q + qbar.
    if (minDauType == 1 && maxDauType == 1 && 
      (motherType == 4 || motherType == 7) ) {
      MEkind = 2;
      if (idMother == 21 || idMother == 22) MEcombi = 1;
      else if (idMother == 23 || idDau1 + idDau2 == 0) {
        MEcombi = 3; 
        dip.MEmix = gammaZmix( event, iMother, dip.iRadiator, dip.iRecoiler );
      }
      else if (idMother == 24) MEcombi = 4;
    }
    // For chi -> chi q qbar, use V/A -> q qbar as first approximation.
    else if (minDauType == 1 && maxDauType == 1 && motherType == 9)
      MEkind = 2;

    // q -> q + V.
    else if (minDauType == 1 && maxDauType == 7 && motherType == 1) 
      MEkind = 3;
      if (idDau1 == 22 || idDau2 == 22) MEcombi = 1;
 
    // Scalar/pseudoscalar -> q + qbar; q -> q + S.
    else if (minDauType == 1 && maxDauType == 1 && motherType == 8) {
      MEkind = 4;
      if (idMother == 25 || idMother == 35 || idMother == 37) MEcombi = 1;
      else if (idMother == 36) MEcombi = 2;
    } 
    else if (minDauType == 1 && maxDauType == 8 && motherType == 1)
      MEkind = 5;
 
    // V -> ~q + ~qbar; ~q -> ~q + V; S -> ~q + ~qbar; ~q -> ~q + S.
    else if (minDauType == 2 && maxDauType == 2 && (motherType == 4 
      || motherType == 7) ) MEkind = 6;
    else if (minDauType == 2 && (maxDauType == 4 || maxDauType == 7) 
      && motherType == 2) MEkind = 7;
    else if (minDauType == 2 && maxDauType == 2 && motherType == 8)
      MEkind = 8;
    else if (minDauType == 2 && maxDauType == 8 && motherType == 2)
      MEkind = 9;
 
    // chi -> q + ~qbar; ~q -> q + chi; q -> ~q + chi.
    else if (minDauType == 1 && maxDauType == 2 && motherType == 9) 
      MEkind = 10;
    else if (minDauType == 1 && maxDauType == 9 && motherType == 2) 
      MEkind = 11;
    else if (minDauType == 2 && maxDauType == 9 && motherType == 1) 
      MEkind = 12;
 
    // ~g -> q + ~qbar; ~q -> q + ~g; q -> ~q + ~g.
    else if (minDauType == 1 && maxDauType == 2 && motherType == 5)
      MEkind = 13;
    else if (minDauType == 1 && maxDauType == 5 && motherType == 2) 
      MEkind = 14;
    else if (minDauType == 2 && maxDauType == 5 && motherType == 1) 
      MEkind = 15;

    // In cases where coloured spin 1 particle involved use spin 0.
    // V_coloured -> q + l.
    else if (minDauType == 1 && maxDauType == 9 && motherType == 3) 
      MEkind = 11; 
    // q -> V_coloured + l;
    else if (minDauType == 3 && maxDauType == 9 && motherType == 1) 
      MEkind = 12;        

    // g (+V, S) -> ~g + ~g (eikonal approximation).
    else if (minDauType == 5 && maxDauType == 5) MEkind = 16;

    // Save ME type and gamma_5 admixture. 
    dip.MEtype = 5 * MEkind + MEcombi; 

  // Now begin processing of charge dipole - still primitive.
  } else if (dip.chgType != 0) {

    // Set defaults for QED case; then possibly done.
    dip.MEorder = true;
    dip.MEsplit = true; 
    if (dip.MEtype >= 0) return;

    // So far only ME corrections for q qbar or l lbar.
    int idDau1 = event[dip.iRadiator].id();
    int idDau2 = event[dip.iMEpartner].id();
    if (abs(idDau1) < 9 && abs(idDau2) < 9 && idDau1 * idDau2 < 0) ;
    else if (abs(idDau1) > 10 && abs(idDau1) < 19 && abs(idDau2) > 10
      && abs(idDau2) < 19 && idDau1 * idDau2 < 0) ;
    else { dip.MEtype = 0; return; }

    // Distinguish charge sum != 0 or = 0; in latter assume vector source.
    dip.MEtype = 101;
    if (idDau1 + idDau2 == 0) dip.MEtype = 102; 
    dip.MEmix = 1.;
  }

}

//--------------------------------------------------------------------------
 
// Find type of particle for ME type: 0 = unknown, 1 = quark, 2 = squark,
// 3 = spare triplet, 4 = gluon, 5 = gluino, 6 = spare octet, 
// 7 = vector boson, 8 = colourless scalar, 9 = colourless spin 1/2.

int TimeShower::findMEparticle( int id, bool isHiddenColour) {

  // find colour and spin of particle.
  int type = 0;
  int colType = abs(particleDataPtr->colType(id)); 
  int spinType = particleDataPtr->spinType(id);

  // For hidden valley particle treat HV colour as normal one.
  // Note: no need to assign gv/gammav since not in ME.
  if (isHiddenColour) {
    colType = 0;
    int idAbs = abs(id);
    if (  (idAbs > 4900000 && idAbs < 4900007)
       || (idAbs > 4900010 && idAbs < 4900017)
       || idAbs == 4900101) colType = 1; 
  } 

  // Find particle type from colour and spin.
  if      (colType == 1 && spinType == 2) type = 1;
  else if (colType == 1 && spinType == 1) type = 2;
  else if (colType == 1)                  type = 3;
  else if (colType == 2 && spinType == 3) type = 4;
  else if (colType == 2 && spinType == 2) type = 5;
  else if (colType == 2)                  type = 6;
  else if (colType == 0 && spinType == 3) type = 7;
  else if (colType == 0 && spinType == 1) type = 8;
  else if (colType == 0 && spinType == 2) type = 9;

  // Done.
  return type;

}  

//--------------------------------------------------------------------------

// Find mixture of V and A in gamma/Z: energy- and flavour-dependent. 

double TimeShower::gammaZmix( Event& event, int iRes, int iDau1, int iDau2) {

  // Try to identify initial flavours; use e+e- as default.
  int idIn1 = -11;
  int idIn2 = 11;
  int iIn1  = (iRes >= 0) ? event[iRes].mother1() : -1;
  int iIn2  = (iRes >= 0) ? event[iRes].mother2() : -1;
  if (iIn1 >=0) idIn1 = event[iIn1].id();
  if (iIn2 >=0) idIn2 = event[iIn1].id();
         
  // In processes f + g/gamma -> f + Z only need find one fermion.
  if (idIn1 == 21 || idIn1 == 22) idIn1 = -idIn2;
  if (idIn2 == 21 || idIn2 == 22) idIn2 = -idIn1;
 
  // Initial flavours and couplings; return if don't make sense.
  if (idIn1 + idIn2 != 0 ) return 0.5;
  int idInAbs = abs(idIn1);
  if (idInAbs == 0 || idInAbs > 18 ) return 0.5; 
  double ei = coupSMPtr->ef(idInAbs);
  double vi = coupSMPtr->vf(idInAbs);
  double ai = coupSMPtr->af(idInAbs);

  // Final flavours and couplings; return if don't make sense.
  if (event[iDau1].id() + event[iDau2].id() != 0) return 0.5;
  int idOutAbs = abs(event[iDau1].id());
  if (idOutAbs == 0 || idOutAbs >18 ) return 0.5; 
  double ef = coupSMPtr->ef(idOutAbs);
  double vf = coupSMPtr->vf(idOutAbs);
  double af = coupSMPtr->af(idOutAbs);

  // Calculate prefactors for interference and resonance part.
  Vec4 psum = event[iDau1].p() + event[iDau2].p();
  double sH = psum.m2Calc();
  double intNorm = 2. * thetaWRat * sH * (sH - mZ*mZ)
    / ( pow2(sH - mZ*mZ) + pow2(sH * gammaZ / mZ) );
  double resNorm = pow2(thetaWRat * sH) 
    / ( pow2(sH - mZ*mZ) + pow2(sH * gammaZ / mZ) );

  // Calculate vector and axial expressions and find mix.
  double vect = ei*ei * ef*ef + ei*vi * intNorm * ef*vf
    + (vi*vi + ai*ai) * resNorm * vf*vf;
  double axiv = (vi*vi + ai*ai) * resNorm * af*af;
  return vect / (vect + axiv);
}

//--------------------------------------------------------------------------

// Set up to calculate QCD ME correction with calcMEcorr.
// Normally for primary particles, but also from g/gamma -> f fbar.
  
double TimeShower::findMEcorr(TimeDipoleEnd* dip, Particle& rad, 
  Particle& partner, Particle& emt, bool cutEdge) {
  
  // Initial values and matrix element kind.
  double wtME    = 1.;
  double wtPS    = 1.; 
  int    MEkind  = dip->MEtype / 5;
  int    MEcombi = dip->MEtype % 5;

  // Construct ME variables.
  Vec4   sum     = rad.p() + partner.p() + emt.p();
  double eCMME   = sum.mCalc();
  double x1      = 2. * (sum * rad.p()) / pow2(eCMME);
  double x2      = 2. * (sum * partner.p()) / pow2(eCMME); 
  double r1      = rad.m() / eCMME;
  double r2      = partner.m() / eCMME; 
  double r3      = 0.;

  // Evaluate kinematics for Hidden Valley with massive gamma_v.
  double gammavCorr = 1.;
  if (dip->colvType != 0 && brokenHVsym) {
    r3              = emt.m() / eCMME;
    double x3Tmp    = 2. - x1 - x2; 
    gammavCorr      = x3Tmp / (x3Tmp - kRad * (x1 + x3Tmp));    
    // For Q_v Qbar_v pair correct kinematics to common average mass.
    if (MEkind == 31) {
      double m2Pair = (rad.p() + partner.p()).m2Calc();
      double m2Avg  = 0.5 * (rad.m2() + partner.m2())  
                    - 0.25 * pow2(rad.m2() - partner.m2()) / m2Pair;
      r1            = sqrt(m2Avg) / eCMME;
      r2            = r1;
      double xShift = 0.5 * (x1 + x2) * (partner.m2() - rad.m2()) / m2Pair;
      x1           += xShift;
      x2           -= xShift;  
    }
  }

  // Derived ME variables, suitably protected.
  double x1minus, x2minus, x3;
  if (cutEdge) {
    x1minus = max(XMARGIN, 1. + r1*r1 - r2*r2 - x1);
    x2minus = max(XMARGIN, 1. + r2*r2 - r1*r1 - x2) ;
    x3      = max(XMARGIN, 2. - x1 - x2);
  } else {
    x1minus = max(XMARGIN*XMARGIN, 1. + r1*r1 - r2*r2 - x1);
    x2minus = max(XMARGIN*XMARGIN, 1. + r2*r2 - r1*r1 - x2) ;
    x3      = max(XMARGIN*XMARGIN, 2. - x1 - x2);
  }

  // Begin processing of QCD dipoles.
  if (dip->colType !=0 || dip->colvType != 0) {

    // Evaluate normal ME, for proper order of particles.
    if (dip->MEorder) wtME = calcMEcorr(MEkind, MEcombi, dip->MEmix, 
      x1, x2, r1, r2, r3, cutEdge);
    else wtME = calcMEcorr(MEkind, MEcombi, dip->MEmix, 
      x2, x1, r2, r1, r3, cutEdge);

    // Split up total ME when two radiating particles.
    if (dip->MEsplit) wtME = wtME * x1minus / x3; 

    // Evaluate shower rate to be compared with.
    wtPS = 2. / ( x3 * x2minus );
    if (dip->MEgluinoRec) wtPS *= 9./4.;
    if (dip->colvType != 0 && brokenHVsym) wtPS *= gammavCorr;
  
  // For generic charge combination currently only massless expression.
  // (Masses included only to respect phase space boundaries.)
  } else if (dip->chgType !=0 && dip->MEtype == 101) {
    double chg1 = particleDataPtr->charge(rad.id());
    double chg2 = particleDataPtr->charge(partner.id());
    wtME = (x1*x1 + x2*x2) * pow2( chg1 * x1minus / x3 
      - chg2 * x2minus / x3 );
    wtPS = 2. * ( chg1*chg1 * x1minus / x3 + chg2*chg2 * x2minus / x3 ); 

  // For flavour neutral system assume vector source and include masses.
  } else if (dip->chgType !=0 && dip->MEtype == 102) {
    wtME = calcMEcorr(2, 1, dip->MEmix, x1, x2, r1, r2, 0., cutEdge) 
      * x1minus / x3;
    wtPS = 2. / ( x3 * x2minus );
  }
  if (wtME > wtPS) infoPtr->errorMsg("Warning in TimeShower::findMEcorr: "
    "ME weight above PS one");
       
  // Return ratio of actual ME to assumed PS rate of emission.
  return wtME / wtPS; 
}

//--------------------------------------------------------------------------

// Matrix elements for gluon (or photon) emission from
// a two-body state; to be used by the parton shower routine.
// Here x_i = 2 E_i/E_cm, r_i = m_i/E_cm and
// 1/sigma_0 d(sigma)/d(x_1)d(x_2) = (alpha-strong/2 pi) * C_F * (this),
// i.e. normalization is such that one recovers the familiar
// (x_1^2 + x_2^2)/((1-x_1)*(1-x_2)) for the massless case.
// Coupling structure:
// kind =  1 : eikonal soft-gluon expression (spin-independent)
//      =  2 : V -> q qbar (V = vector/axial vector colour singlet)
//      =  3 : q -> q V
//      =  4 : S -> q qbar (S = scalar/pseudoscalar colour singlet)
//      =  5 : q -> q S
//      =  6 : V -> ~q ~qbar (~q = squark)
//      =  7 : ~q -> ~q V
//      =  8 : S -> ~q ~qbar
//      =  9 : ~q -> ~q S
//      = 10 : chi -> q ~qbar (chi = neutralino/chargino)
//      = 11 : ~q -> q chi
//      = 12 : q -> ~q chi
//      = 13 : ~g -> q ~qbar
//      = 14 : ~q -> q ~g
//      = 15 : q -> ~q ~g
//      = 16 : (9/4)*(eikonal) for gg -> ~g ~g
//      = 30 : Dv -> d qv     (Dv= hidden valley fermion, qv= valley scalar)
//      = 31 : S  -> Dv Dvbar (S=scalar color singlet)
// Note that the order of the decay products is important.
// combi = 1 : pure non-gamma5, i.e. vector/scalar/...
//       = 2 : pure gamma5, i.e. axial vector/pseudoscalar/....
//       = 3 : mixture mix*(combi=1) + (1-mix)*(combi=2)
//       = 4 : mixture (combi=1) +- (combi=2)

double TimeShower::calcMEcorr( int kind, int combiIn, double mixIn, 
  double x1, double x2, double r1, double r2, double r3, bool cutEdge) {

  // Frequent variable combinations.
  double x3     = 2. - x1 - x2;
  double x1s    = x1 * x1;
  double x2s    = x2 * x2;
  double x3s    = x3 * x3;
  double x1c    = x1 * x1s;
  double x2c    = x2 * x2s;
  double x3c    = x3 * x3s;
  double r1s    = r1 * r1;
  double r2s    = r2 * r2;
  double r1c    = r1 * r1s;
  double r2c    = r2 * r2s;
  double r1q    = r1s * r1s;
  double r2q    = r2s * r2s;
  double prop1  = 1. + r1s - r2s - x1; 
  double prop2  = 1. + r2s - r1s - x2;
  double prop1s = prop1 * prop1;
  double prop2s = prop2 * prop2;
  double prop12 = prop1 * prop2;
  double prop13 = prop1 * x3;
  double prop23 = prop2 * x3;

  // Special case: Hidden-Valley massive photon. 
  double r3s    = r3 * r3;
  double prop3  = r3s - x3;
  double prop3s = prop3 * prop3;
  if (kind == 30) prop13 = prop1 * prop3;

  // Check input values. Return zero outside allowed phase space.
  if (cutEdge) {
    if (x1 - 2.*r1 < XMARGIN || prop1 < XMARGIN) return 0.;
    if (x2 - 2.*r2 < XMARGIN || prop2 < XMARGIN) return 0.;
    // Limits not worked out for r3 > 0.
    if (kind != 30 && kind != 31) {
      if (x1 + x2 - 1. - pow2(r1+r2) < XMARGIN) return 0.;
      // Note: equivalent rewritten form 4. * ( (1. - x1) * (1. - x2) 
      // * (1. - r1s - r2s - x3) + r1s * (1. - x2s - x3) + r2s 
      // * (1. - x1s - x3) - pow2(r1s - r2s) ) gives about same result.
      if ( (x1s - 4.*r1s) * (x2s - 4.*r2s) 
        - pow2( 2. * (1. - x1 - x2 + r1s + r2s) + x1*x2 ) 
        < XMARGIN * (XMARGINCOMB + r1 + r2) ) return 0.;
    }
  }

  // Initial values; phase space.
  int combi   = max(1, min(4, combiIn) ); 
  double mix  = max(0., min(1., mixIn) );
  bool isSet1 = false;
  bool isSet2 = false;
  bool isSet4 = false;
  double ps = sqrtpos( pow2(1. - r1*r1 - r2*r2) - pow2(2. * r1 * r2) );
  double rLO = 0., rFO = 0., rLO1 = 0., rFO1 = 0., rLO2 = 0., 
    rFO2 = 0., rLO4 = 0., rFO4 = 0.;
  double offset = 0;
 
  // Select which kind of ME to use.
  switch (kind) {

    // case 1 is equal to default, i.e. eikonal expression.

    // V -> q qbar (V = gamma*/Z0/W+-/...).
    case 2:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(2.-r1s-r1q+6.*r1*r2-r2s+2.*r1s*r2s-r2q)/2.;
        rFO1 = -(3.+6.*r1s+r1q-6.*r1*r2+6.*r1c*r2-2.*r2s-6.*r1s*r2s
        +6.*r1*r2c+r2q-3.*x1+6.*r1*r2*x1+2.*r2s*x1+x1s-2.*r1s*x1s
        +3.*r1s*x3+6.*r1*r2*x3-r2s*x3-2.*x1*x3-5.*r1s*x1*x3
        +r2s*x1*x3+x1s*x3-3.*x3s-3.*r1s*x3s+r2s*x3s
        +2.*x1*x3s+x3c-x2)
        /prop2s
        -2.*(-3.+r1s-6.*r1*r2+6.*r1c*r2+3.*r2s-4.*r1s*r2s
        +6.*r1*r2c+2.*x1+3.*r1s*x1+r2s*x1-x1s-r1s*x1s
        -r2s*x1s+4.*x3+2.*r1s*x3+3.*r1*r2*x3-r2s*x3-3.*x1*x3
        -2.*r1s*x1*x3+x1s*x3-x3s-r1s*x3s+r1*r2*x3s+x1*x3s)
        /prop12
        -(-1.+2.*r1s+r1q+6.*r1*r2+6.*r1c*r2-2.*r2s-6.*r1s*r2s
        +6.*r1*r2c+r2q-x1-2.*r1s*x1-6.*r1*r2*x1+8.*r2s*x1+x1s
        -2.*r2s*x1s-r1s*x3+r2s*x3-r1s*x1*x3+r2s*x1*x3+x1s*x3+x2)
        /prop1s;
        rFO1 = rFO1/2.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(2.-r1s-r1q-6.*r1*r2-r2s+2.*r1s*r2s-r2q)/2.;
        rFO2 = -(3.+6.*r1s+r1q+6.*r1*r2-6.*r1c*r2-2.*r2s-6.*r1s*r2s    
        -6.*r1*r2c+r2q-3.*x1-6.*r1*r2*x1+2.*r2s*x1+x1s-2.*r1s*x1s
        +3.*r1s*x3-6.*r1*r2*x3-r2s*x3-2.*x1*x3-5.*r1s*x1*x3
        +r2s*x1*x3+x1s*x3-3.*x3s-3.*r1s*x3s+r2s*x3s+2.*x1*x3s+x3c-x2)
        /prop2s
        -2.*(-3+r1s+6.*r1*r2-6.*r1c*r2+3.*r2s-4.*r1s*r2s-6.*r1*r2c
        +2.*x1+3.*r1s*x1+r2s*x1-x1s-r1s*x1s-r2s*x1s+4.*x3+2.*r1s*x3
        -3.*r1*r2*x3-r2s*x3-3.*x1*x3-2.*r1s*x1*x3+x1s*x3-x3s-r1s*x3s
        -r1*r2*x3s+x1*x3s)
        /prop12
        -(-1.+2.*r1s+r1q-6.*r1*r2-6.*r1c*r2-2.*r2s-6.*r1s*r2s
        -6.*r1*r2c+r2q-x1-2.*r1s*x1+6.*r1*r2*x1+8.*r2s*x1+x1s
        -2.*r2s*x1s-r1s*x3+r2s*x3-r1s*x1*x3+r2s*x1*x3+x1s*x3+x2)
        /prop1s;
        rFO2 = rFO2/2.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(2.-r1s-r1q-r2s+2.*r1s*r2s-r2q)/2.;
        rFO4 = (1.-r1q+6.*r1s*r2s-r2q+x1+3.*r1s*x1-9.*r2s*x1-3.*x1s
        -r1s*x1s+3.*r2s*x1s+x1c-x2-r1s*x2+r2s*x2-r1s*x1*x2+r2s*x1*x2
        +x1s*x2)
        /prop1s 
        -2.*(1.+r1s+r2s-4.*r1s*r2s+r1s*x1+2.*r2s*x1-x1s-r2s*x1s
        +2.*r1s*x2+r2s*x2-3.*x1*x2+x1s*x2-x2s-r1s*x2s+x1*x2s)
        /prop12
        +(1.-r1q+6.*r1s*r2s-r2q-x1+r1s*x1-r2s*x1+x2-9.*r1s*x2
        +3.*r2s*x2+r1s*x1*x2-r2s*x1*x2-3.*x2s+3.*r1s*x2s-r2s*x2s
        +x1*x2s+x2c)
        /prop2s;
        rFO4 = rFO4/2.;
        isSet4 = true;
      }
      break; 
 
    // q -> q V.
    case 3:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-2.*r1s+r1q+r2s-6.*r1*r2s+r1s*r2s-2.*r2q);
        rFO1 = -2.*(-1.+r1-2.*r1s+2.*r1c-r1q+pow5(r1)-r2s+r1*r2s
        -5.*r1s*r2s+r1c*r2s-2.*r1*r2q+2.*x1-2.*r1*x1+2.*r1s*x1
        -2.*r1c*x1+2.*r2s*x1+5.*r1*r2s*x1+r1s*r2s*x1+2.*r2q*x1
        -x1s+r1*x1s-r2s*x1s+3.*x2+4.*r1s*x2+r1q*x2+2.*r2s*x2
        +2.*r1s*r2s*x2-4.*x1*x2-2.*r1s*x1*x2-r2s*x1*x2+x1s*x2
        -2.*x2s-2.*r1s*x2s+x1*x2s)
        /prop23
        +(2.*r2s+6.*r1*r2s-6.*r1s*r2s+6.*r1c*r2s+2.*r2q+6.*r1*r2q
        -r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2-3.*r2s*x2-6.*r1*r2s*x2
        +9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1+6.*r1*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s-
        2.*r1s*x1s+x1c+7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2+6.*r1*r2s*x2
        +r1s*r2s*x2-2.*r2q*x2-9.*x1*x2-3.*r1s*x1*x2+2.*r2s*x1*x2
        +2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s+x1*x2s)
        /x3s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-2.*r1s+r1q+r2s+6.*r1*r2s+r1s*r2s-2.*r2q);
        rFO2 = 2*(1.+r1+2.*r1s+2.*r1c+r1q+pow5(r1)+r2s+r1*r2s
        +5.*r1s*r2s+r1c*r2s-2.*r1*r2q-2.*x1-2.*r1*x1-2.*r1s*x1
        -2.*r1c*x1-2.*r2s*x1+5.*r1*r2s*x1-r1s*r2s*x1-2.*r2q*x1+x1s
        +r1*x1s+r2s*x1s-3.*x2-4.*r1s*x2-r1q*x2-2.*r2s*x2
        -2.*r1s*r2s*x2+4.*x1*x2+2.*r1s*x1*x2+r2s*x1*x2-x1s*x2
        +2.*x2s+2.*r1s*x2s-x1*x2s)
        /prop23
        +(2.*r2s-6.*r1*r2s-6.*r1s*r2s-6.*r1c*r2s+2.*r2q-6.*r1*r2q
        -r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2-3.*r2s*x2+6.*r1*r2s*x2
        +9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1-6.*r1*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s
        -2.*r1s*x1s+x1c+7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2-6.*r1*r2s*x2
        +r1s*r2s*x2-2.*r2q*x2-9.*x1*x2-3.*r1s*x1*x2+2.*r2s*x1*x2
        +2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s+x1*x2s)
        /x3s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-2.*r1s+r1q+r2s+r1s*r2s-2.*r2q);
        rFO4 = 2*(1.+2.*r1s+r1q+r2s+5.*r1s*r2s-2.*x1-2.*r1s*x1
        -2.*r2s*x1-r1s*r2s*x1-2.*r2q*x1+x1s+r2s*x1s-3.*x2-4.*r1s*x2
        -r1q*x2-2.*r2s*x2-2.*r1s*r2s*x2+4.*x1*x2+2.*r1s*x1*x2+r2s*x1*x2
        -x1s*x2+2.*x2s+2.*r1s*x2s-x1*x2s)
        /prop23
        +(2.*r2s-6.*r1s*r2s+2.*r2q-r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2
        -3.*r2s*x2+9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s-2.*r1s*x1s+x1c
        +7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2+r1s*r2s*x2-2.*r2q*x2-9.*x1*x2
        -3.*r1s*x1*x2+2.*r2s*x1*x2+2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s
        +x1*x2s)
        /x3s;
        isSet4 = true;
      }
      break; 
 
    // S -> q qbar    (S = h0/H0/A0/H+-/...).
    case 4:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s-r2s-2.*r1*r2);
        rFO1 = -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -2.*(r1s+r1q-2.*r1c*r2+r2s-6.*r1s*r2s-2.*r1*r2c+r2q-r1s*x1
        +r1*r2*x1+2.*r2s*x1+2.*r1s*x2+r1*r2*x2-r2s*x2-x1*x2)
        /prop12
        -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1-r1s*x1
        +r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s-r2s+2.*r1*r2);
        rFO2 = -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1-2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +2.*(-r1s-r1q-2.*r1c*r2-r2s+6.*r1s*r2s-2.*r1*r2c-r2q+r1s*x1
        +r1*r2*x1-2.*r2s*x1-2.*r1s*x2+r1*r2*x2+r2s*x2+x1*x2)
        /prop12;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+3.*r2s*x1+x2
        +r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -2.*(r1s+r1q+r2s-6.*r1s*r2s+r2q-r1s*x1
        +2.*r2s*x1+2.*r1s*x2-r2s*x2-x1*x2)
        /prop12
        -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1
        +x2+3.*r1s*x2-r2s*x2-x1*x2)
        /prop2s;
        isSet4 = true;
      }
      break; 
 
    // q -> q S.
    case 5:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = (4.-4.*r1s+4.*r2s-3.*x1-2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        -2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.-r1-5.*r1s-r1c+3.*r2s+r1*r2s-2.*x1-r1*x1
        +r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.+r1s-r2s-2.*r1);
        rFO2 = (4.-4.*r1s+4.*r2s-3.*x1+2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        +2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.+r1-5.*r1s+r1c+3.*r2s-r1*r2s-2.*x1+r1*x1
        +r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = (4.-4.*r1s+4.*r2s-3.*x1+r1s*x1-r2s*x1-5.*x2+r1s*x2
        -r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.-5.*r1s+3.*r2s-2.*x1+r1s*x1-4.*x2+2.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop23
        +(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet4 = true;
      }
      break; 
 
    // V -> ~q ~qbar  (~q = squark).
    case 6:
      rLO1 = ps*(1.-2.*r1s+r1q-2.*r2s-2.*r1s*r2s+r2q);
      rFO1 = 2.*3.+(1.+r1s+r2s-x1)*(4.*r1s-x1s)
      /prop1s
      +2.*(-1.-3.*r1s-r2s+x1+x1s*0.5+x2-x1*x2*0.5)
      /prop1
      +(1.+r1s+r2s-x2)*(4.*r2s-x2s)
      /prop2s
      +2.*(-1.-r1s-3.*r2s+x1+x2-x1*x2*0.5+x2s*0.5)
      /prop2
      -(-4.*r1s-4.*r1q-4.*r2s-8.*r1s*r2s-4.*r2q+2.*x1+6.*r1s*x1
      +6.*r2s*x1-2.*x1s+2.*x2+6.*r1s*x2+6.*r2s*x2-4.*x1*x2
      -2.*r1s*x1*x2-2.*r2s*x1*x2+x1s*x2-2.*x2s+x1*x2s)
      /prop12;
      isSet1 = true;
      break; 
 
    // ~q -> ~q V.
    case 7:
      rLO1 = ps*(1.-2.*r1s+r1q-2.*r2s-2.*r1s*r2s+r2q);
      rFO1 = 16.*r2s-8.*(4.*r2s+2.*r2s*x1+x2+r1s*x2+r2s*x2-x1*x2
      -2.*x2s)
      /(3.*prop2)
      +8.*(1.+r1s+r2s-x2)*(4.*r2s-x2s)
      /(3.*prop2s)
      +8.*(x1+x2)*(-1.-2.*r1s-r1q-2.*r2s+2.*r1s*r2s-r2q+2.*x1
      +2.*r1s*x1+2.*r2s*x1-x1s+2.*x2+2.*r1s*x2+2.*r2s*x2-2.*x1*x2-x2s)
      /(3.*x3s)
      +8.*(-1.-r1s+r2s-x1)*(2.*r2s*x1+x2+r1s*x2+r2s*x2-x1*x2-x2s)
      /(3.*prop2*x3)
      -8.*(1.+2.*r1s+r1q+2.*r2s-2.*r1s*r2s+r2q-2.*x1-2.*r1s*x1
      -4.*r2s*x1+x1s-3.*x2-3.*r1s*x2-3.*r2s*x2+3.*x1*x2+2.*x2s)
      /(3.*x3);
      rFO1 = 3.*rFO1/8.;
      isSet1 = true;
      break; 
 
    // S -> ~q ~qbar.
    case 8:
      rLO1 = ps;
      rFO1 = (-1.-2.*r1s-r1q-2.*r2s+2.*r1s*r2s-r2q+2.*x1+2.*r1s*x1
      +2.*r2s*x1-x1s-r2s*x1s+2.*x2+2.*r1s*x2+2.*r2s*x2-3.*x1*x2
      -r1s*x1*x2-r2s*x1*x2+x1s*x2-x2s-r1s*x2s+x1*x2s)
      /(prop1s*prop2s);
      rFO1 = 2.*rFO1;
      isSet1 = true;
      break; 
 
    // ~q -> ~q S.
    case 9:
      rLO1 = ps;
      rFO1 = (-1.-r1s-r2s+x2)
      /prop2s
      +(1.+r1s-r2s+x1)
      /prop23
      -(x1+x2)
      /x3s;
      isSet1 = true;
      break; 
 
    // chi -> q ~qbar   (chi = neutralino/chargino).
    case 10:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = (2.*r1+x1)*(-1.-r1s-r2s+x1)
        /prop1s
        +2.*(-1.-r1s-2.*r1c-r2s-2.*r1*r2s+3.*x1*0.5+r1*x1
        -r1s*x1*0.5-r2s*x1*0.5+x2+r1*x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-2.*r1+r1s-r2s);
        rFO2 = (2.*r1-x1)*(1.+r1s+r2s-x1)
        /prop1s
        +2.*(-1.-r1s+2.*r1c-r2s+2.*r1*r2s+3.*x1*0.5-r1*x1
        -r1s*x1*0.5-r2s*x1*0.5+x2-r1*x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)/
        prop2s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = x1*(-1.-r1s-r2s+x1)
        /prop1s
        +2.*(-1.-r1s-r2s+3.*x1*0.5-r1s*x1*0.5-r2s*x1*0.5
        +x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet4 = true;
      }
      break; 
 
    // ~q -> q chi.
    case 11:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-pow2(r1+r2));
        rFO1 = (1.+r1s+2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q-2.*r1*r2-2.*r1c*r2+2.*r1*r2c+r2q+x1+r1s*x1
        -2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-pow2(r1-r2));
        rFO2 = (1.+r1s-2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q+2.*r1*r2+2.*r1c*r2-2.*r1*r2c+r2q+x1+r1s*x1
        +2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = (1.+r1s+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1+x2
        +3.*r1s*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q+r2q+x1+r1s*x1-3.*r2s*x1
        +2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet4 = true;
      }
      break; 
 
    // q -> ~q chi.
    case 12:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s+r2s+2.*r2);
        rFO1 = (2.*r2+x2)*(-1.-r1s-r2s+x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1-2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2-2.*r2*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s+r2+r1s*r2-r2s-r2c+x1+r2*x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s+r2s-2.*r2);
        rFO2 = (2.*r2-x2)*(1.+r1s+r2s-x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2+2.*r2*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s-r2-r1s*r2-r2s+r2c+x1-r2*x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s+r2s);
        rFO4 = x2*(-1.-r1s-r2s+x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+r2s*x1+x1s
        -3.*x2-r1s*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s-r2s+x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet4 = true;
      }
      break; 
 
    // ~g -> q ~qbar.
    case 13:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = 4.*(2.*r1+x1)*(-1.-r1s-r2s+x1)
        /(3.*prop1s)
        -(-1.-r1s-2.*r1c-r2s-2.*r1*r2s+3.*x1*0.5+r1*x1-r1s*x1*0.5
        -r2s*x1*0.5+x2+r1*x2+r1s*x2-x1*x2*0.5)
        /(3.*prop12)
        +3.*(-1.+r1-r1s-r1c-r2s+r1*r2s+2.*x1+r2s*x1-x1s*0.5+x2+r1*x2
        +r1s*x2-x1*x2*0.5)
        /prop13
        +3.*(4.-4.*r1s+4.*r2s-3.*x1-2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        -2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -3.*(3.-r1-5.*r1s-r1c+3.*r2s+r1*r2s-2.*x1-r1*x1+r1s*x1
        -4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +4.*(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1-r2s*x1
        -3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /(3.*prop2s);
        rFO1 = 3.*rFO1/4.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.+r1s-r2s-2.*r1);
        rFO2 = 4.*(2.*r1-x1)*(1.+r1s+r2s-x1)
        /(3.*prop1s)
        +3.*(-1.-r1-r1s+r1c-r2s-r1*r2s+2.*x1+r2s*x1-x1s*0.5
        +x2-r1*x2+r1s*x2-x1*x2*0.5)
        /prop13
        +(2.+2.*r1s-4.*r1c+2.*r2s-4.*r1*r2s-3.*x1+2.*r1*x1
        +r1s*x1+r2s*x1-2.*x2+2.*r1*x2-2.*r1s*x2+x1*x2)
        /(6.*prop12)
        +3.*(4.-4.*r1s+4.*r2s-3.*x1+2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        +2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -3.*(3.+r1-5.*r1s+r1c+3.*r2s-r1*r2s-2.*x1+r1*x1+r1s*x1-4.*x2
        +2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +4.*(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1-r2s*x1
        -3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /(3.*prop2s);
        rFO2 = 3.*rFO2/4.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = 8.*x1*(-1.-r1s-r2s+x1)
        /(3.*prop1s)
        +6.*(-1-r1s-r2s+2.*x1+r2s*x1-x1s*0.5+x2+r1s*x2-x1*x2*0.5)
        /prop13
        +(2.+2.*r1s+2.*r2s-3.*x1+r1s*x1+r2s*x1-2.*x2-2.*r1s*x2+x1*x2)
        /(3.*prop12)
        +6.*(4.-4.*r1s+4.*r2s-3.*x1+r1s*x1-r2s*x1-5.*x2+r1s*x2-r2s*x2
        +x1*x2+x2s)
        /x3s
        -6.*(3.-5.*r1s+3.*r2s-2.*x1+r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +8.*(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2-r2s*x2
        +x1*x2+x2s)
        /(3.*prop2s);
        rFO4 = 3.*rFO4/8.;
        isSet4 = true;
      }
      break; 
 
    // ~q -> q ~g.
    case 14:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s-r2s-2.*r1*r2);
        rFO1 = 64.*(1.+r1s+2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -16.*(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q
        +x1-r1s*x1+2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -16.*(r1s+r1q-2.*r1c*r2+r2s-6.*r1s*r2s-2.*r1*r2c+r2q-r1s*x1
        +r1*r2*x1+2.*r2s*x1+2.*r1s*x2+r1*r2*x2-r2s*x2-x1*x2)
        /prop12
        -64.*(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /(9.*prop2s)
        +8.*(-1.+r1q-2.*r1*r2+2.*r1c*r2-2.*r2s-2.*r1*r2c-r2q-2.*r1s*x1
        +2.*r2s*x1+x1s+x2-3.*r1s*x2-2.*r1*r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-2.*r1s-r1q-2.*r1*r2-2.*r1c*r2+2.*r1*r2c+r2q+x1+r1s*x1
        -2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO1 = 9.*rFO1/64.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s-r2s+2.*r1*r2);
        rFO2 = 64.*(1.+r1s-2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -16.*(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1-2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -64.*(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /(9.*prop2s)
        +16.*(-r1s-r1q-2.*r1c*r2-r2s+6.*r1s*r2s-2.*r1*r2c-r2q+r1s*x1
        +r1*r2*x1-2.*r2s*x1-2.*r1s*x2+r1*r2*x2+r2s*x2+x1*x2)
        /prop12
        +8.*(-1.+r1q+2.*r1*r2-2.*r1c*r2-2.*r2s+2.*r1*r2c-r2q-2.*r1s*x1
        +2.*r2s*x1+x1s+x2-3.*r1s*x2+2.*r1*r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-2.*r1s-r1q+2.*r1*r2+2.*r1c*r2-2.*r1*r2c+r2q+x1+r1s*x1+
        2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO2 = 9.*rFO2/64.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = 128.*(1.+r1s+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -32*(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+3.*r2s*x1+x2
        +r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -32.*(r1s+r1q+r2s-6.*r1s*r2s+r2q-r1s*x1+2.*r2s*x1+2.*r1s*x2
        -r2s*x2-x1*x2)
        /prop12
        -128.*(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1+x2+3.*r1s*x2
        -r2s*x2-x1*x2)
        /(9.*prop2s)
        +16.*(-1.+r1q-2.*r2s-r2q-2.*r1s*x1+2.*r2s*x1+x1s
        +x2-3.*r1s*x2+r2s*x2+x1*x2)
        /prop13
        -16.*(-1.-2.*r1s-r1q+r2q+x1+r1s*x1-3.*r2s*x1
        +2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO4 = 9.*rFO4/128.;
        isSet4 = true;
      }
      break; 
 
    // q -> ~q ~g.
    case 15:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s+r2s+2.*r2);
        rFO1 = 32*(2.*r2+x2)*(-1.-r1s-r2s+x2)
        /(9.*prop2s)
        +8.*(-1.-r1s-2.*r1s*r2-r2s-2.*r2c+x1+r2*x1+r2s*x1
        +3.*x2*0.5-r1s*x2*0.5+r2*x2-r2s*x2*0.5-x1*x2*0.5)
        /prop12
        +8.*(2.+2.*r1s-2.*r2-2.*r1s*r2-6.*r2s-2.*r2c-3.*x1-r1s*x1
        +2.*r2*x1+3.*r2s*x1+x1s-x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        +32.*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1-2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2-2.*r2*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        -8.*(3.+3.*r1s-r2+r1s*r2-5.*r2s-r2c-4.*x1-r1s*x1
        +2.*r2s*x1+x1s-2.*x2-r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-r1s+r2+r1s*r2-r2s-r2c+x1+r2*x1+r2s*x1+2.*x2+r1s*x2
        -x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO1 = 9.*rFO1/32.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s+r2s-2.*r2);
        rFO2 = 32*(2.*r2-x2)*(1.+r1s+r2s-x2)
        /(9.*prop2s)
        +8.*(-1.-r1s+2.*r1s*r2-r2s+2.*r2c+x1-r2*x1+r2s*x1
        +3.*x2*0.5-r1s*x2*0.5-r2*x2-r2s*x2*0.5-x1*x2*0.5)
        /prop12
        +8.*(2.+2.*r1s+2.*r2+2.*r1s*r2-6.*r2s+2.*r2c-3.*x1-r1s*x1
        -2.*r2*x1+3.*r2s*x1+x1s-x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        -8.*(3.+3.*r1s+r2-r1s*r2-5.*r2s+r2c-4.*x1-r1s*x1+2.*r2s*x1+x1s
        -2.*x2+r2*x2+r2s*x2+x1*x2)
        /prop13
        +32*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+2.*r2*x1+r2s*x1
        +x1s-3.*x2-r1s*x2+2.*r2*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        -8.*(-1.-r1s-r2-r1s*r2-r2s+r2c+x1-r2*x1+r2s*x1+2.*x2+r1s*x2
        -x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO2 = 9.*rFO2/32.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s+r2s);
        rFO4 = 64.*x2*(-1.-r1s-r2s+x2)
        /(9.*prop2s)
        +16.*(-1.-r1s-r2s+x1+r2s*x1+3.*x2*0.5-r1s*x2*0.5
        -r2s*x2*0.5-x1*x2*0.5)
        /prop12
        -16.*(3.+3.*r1s-5.*r2s-4.*x1-r1s*x1+2.*r2s*x1+x1s-2.*x2+r2s*x2
        +x1*x2)
        /prop13
        +64.*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+r2s*x1+x1s-3.*x2
        -r1s*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        +16.*(2.+2.*r1s-6.*r2s-3.*x1-r1s*x1+3.*r2s*x1+x1s
        -x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        -16.*(-1.-r1s-r2s+x1+r2s*x1+2.*x2+r1s*x2-x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO4 = 9.*rFO4/64.;
        isSet4 = true;
      }
      break; 
 
    // g -> ~g ~g. Use (9/4)*eikonal. May be changed in the future.
    case 16:
      rLO = ps;
      if      (combi == 2) offset = x3s;
      else if (combi == 3) offset = mix * x3s;
      else if (combi == 4) offset = 0.5 * x3s;
      rFO = ps * 4.5 * ( (x1+x2-1.+offset-r1s-r2s)/prop12 
      - r1s/prop2s - r2s/prop1s );
      break; 

    // Dv -> qv d.
    case 30:
      rLO = ps*(1.-r1s+r2s+2.*r2);
      rFO = ( 0.5*r3s + 2.*r1q + 0.5*r2s*r3s + r2*r3s - 2.*r1s 
             - 0.5*r1s*r3s - 2.*r1s*r2s - 4.*r1s*r2 ) / prop2s
          + ( -2. + 2.*r2q + 2.*r1q + 2.*r2s*r3s - 4.*r2 + 2.*r2*r3s 
             + 4.*r2*r2s - 4.*r1s*r2s - 4.*r1s*r2 ) /prop23
          + ( -2. - 0.5*r3s - 2.*r2s - 4.*r2 + 2.*r1s ) / prop2
          + ( -2. - r3s - 2.*r2s - r2s*r3s - 4.*r2 - 2.*r2*r3s 
             + 2.*r1s + r1s*r3s ) / prop3s
          + ( -1. - r3s - r2s - 4.*r2 + r1s - x2 ) / prop3
          + 1.;
      break;

    // S -> Dv Dvbar
    case 31:
      rLO = ps*(1.-4.*r1s);
      rFO = (r3s + 2.*r1s) * (-1. + 4.*r1s) * (1./prop1s + 1./prop2s)
          + (-1. + 8.*r1s - x2) / prop1
          + (-1. + 8.*r1s - x1) / prop2          
          + 2. * (1. - 6.*r1s + 8.*r1q + 4.*r3s*r1s) / prop12
          + 2.;
      break;

    // Eikonal expression for kind == 1; also acts as default.
    default:
      rLO = ps;
      if      (combi == 2) offset = x3s;
      else if (combi == 3) offset = mix * x3s;
      else if (combi == 4) offset = 0.5 * x3s;
      rFO = ps * 2. * ( (x1+x2-1.+offset-r1s-r2s)/prop12 
      - r1s/prop2s - r2s/prop1s );
      break;

  // End of ME cases. 
  }

  // Find relevant leading and first order expressions.
  if      (combi == 1 && isSet1) {
    rLO = rLO1; 
    rFO = rFO1; }     
  else if (combi == 2 && isSet2) {
    rLO = rLO2; 
    rFO = rFO2; }     
  else if (combi == 3 && isSet1 && isSet2) {
    rLO = mix * rLO1 + (1.-mix) * rLO2; 
    rFO = mix * rFO1 + (1.-mix) * rFO2; }
  else if (isSet4) {
    rLO = rLO4; 
    rFO = rFO4; }     
  else if (combi == 4 && isSet1 && isSet2) {
    rLO = 0.5 * (rLO1 + rLO2);
    rFO = 0.5 * (rFO1 + rFO2); }
  else if (isSet1) {
    rLO = rLO1; 
    rFO = rFO1; } 

  // Return ratio of first to leading order cross section.     
  return rFO / rLO;
}  

//--------------------------------------------------------------------------

// Find coefficient of azimuthal asymmetry from gluon polarization.

void TimeShower::findAsymPol( Event& event, TimeDipoleEnd* dip) {

  // Default is no asymmetry. Only gluons are studied.
  dip->asymPol = 0.;
  dip->iAunt = 0;
  int iRad = dip->iRadiator;
  if (!doPhiPolAsym || event[iRad].id() != 21) return;

  // Trace grandmother via possibly intermediate recoil copies.
  int iMother = event.iTopCopy(iRad);
  int iGrandM = event[iMother].mother1();

  // If grandmother in initial state of hard scattering,
  // then only keep gg and qq initial states.
  int statusGrandM = event[iGrandM].status();
  bool isHardProc  = (statusGrandM == -21 || statusGrandM == -31);  
  if (isHardProc) {
    if (event[iGrandM + 1].status() != statusGrandM) return;
    if (event[iGrandM].isGluon() && event[iGrandM + 1].isGluon());
    else if (event[iGrandM].isQuark() && event[iGrandM + 1].isQuark());
    else return;
  }

  // Set aunt by history or, for hard scattering, by colour flow.
  if (isHardProc) dip->iAunt = dip->iRecoiler;
  else dip->iAunt = (event[iGrandM].daughter1() == iMother) 
    ? event[iGrandM].daughter2() : event[iGrandM].daughter1();

  // Coefficient from gluon production (approximate z by energy).
  // For hard process arbitrarily put z = 1/2.
  double zProd = (isHardProc) ? 0.5 : event[iRad].e() 
    / (event[iRad].e() + event[dip->iAunt].e());
  if (event[iGrandM].isGluon()) dip->asymPol = pow2( (1. - zProd) 
    / (1. - zProd * (1. - zProd) ) );
  else dip->asymPol = 2. * (1. - zProd) / (1. + pow2(1. - zProd) );

  // Coefficients from gluon decay.
  if (dip->flavour == 21) dip->asymPol *= pow2( (1. - dip->z) 
    / (1. - dip->z * (1. - dip->z) ) );
  else  dip->asymPol *= -2. * dip->z *( 1. - dip->z ) 
    / (1. - 2. * dip->z * (1. - dip->z) );

}

//--------------------------------------------------------------------------

// Print the list of dipoles.

void TimeShower::list(ostream& os) const {

  // Header.
  os << "\n --------  PYTHIA TimeShower Dipole Listing  ----------------"
     << "--------------------------------------------- \n \n    i    rad"
     << "    rec       pTmax  col  chg  gam  oni   hv  isr  sys sysR typ"
     << "e  MErec     mix  ord  spl  ~gR \n" << fixed << setprecision(3);
  
  // Loop over dipole list and print it.
  for (int i = 0; i < int(dipEnd.size()); ++i) 
  os << setw(5) << i                     << setw(7) << dipEnd[i].iRadiator 
     << setw(7) << dipEnd[i].iRecoiler   << setw(12) << dipEnd[i].pTmax 
     << setw(5) << dipEnd[i].colType     << setw(5) << dipEnd[i].chgType
     << setw(5) << dipEnd[i].gamType     << setw(5) << dipEnd[i].isOctetOnium 
     << setw(5) << dipEnd[i].isHiddenValley << setw(5) << dipEnd[i].isrType 
     << setw(5) << dipEnd[i].system      << setw(5) << dipEnd[i].systemRec
     << setw(5) << dipEnd[i].MEtype      << setw(7) << dipEnd[i].iMEpartner
     << setw(8) << dipEnd[i].MEmix       << setw(5) << dipEnd[i].MEorder
     << setw(5) << dipEnd[i].MEsplit     << setw(5) << dipEnd[i].MEgluinoRec 
     << "\n";
 
  // Done.
  os << "\n --------  End PYTHIA TimeShower Dipole Listing  ------------"
     << "---------------------------------------------" << endl;
  
}

//==========================================================================

} // end namespace Pythia8

