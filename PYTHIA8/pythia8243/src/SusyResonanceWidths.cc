// SusyResonanceWidths.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand
// Authors: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// SUSY Resonance widths classes.

#include "Pythia8/SusyResonanceWidths.h"
#include "Pythia8/SusyWidthFunctions.h"
#include "Pythia8/SusyCouplings.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaComplex.h"

namespace Pythia8 {

//==========================================================================

// The SUSYResonanceWidths Class
// Derived class for SUSY resonances

const bool SUSYResonanceWidths::DBSUSY = false;

//--------------------------------------------------------------------------

bool SUSYResonanceWidths::initBSM(){

  if (couplingsPtr->isSUSY) {
    coupSUSYPtr     = (CoupSUSY *) couplingsPtr;
    return true;
  }

  return false;
}

//--------------------------------------------------------------------------

bool SUSYResonanceWidths::allowCalc(){

  // Check if decay calculations at all possible
  if ( !couplingsPtr->isSUSY ) return false;
  if ( (idRes == 45 || idRes == 46 || idRes == 1000045)
       && !coupSUSYPtr->isNMSSM ) return false;

  if (settingsPtr->flag("SLHA:useDecayTable") ) {

    // Next check if decay table was read in via SLHA and takes precedence
    for ( int iDec = 0; iDec < int((coupSUSYPtr->slhaPtr)->decays.size());
          ++iDec)
      if ( (coupSUSYPtr->slhaPtr)->decays[iDec].getId() == abs(idRes) ) {
        if (DBSUSY) cout<<"Using external decay table for:"<<idRes<<endl;
        return false;
      }

  }

  // Else we should do the calculation; set available channels
  bool done = getChannels(idRes);
  stringstream idStream;
  idStream << "ID = " << idRes ;
  if (!done)  infoPtr->errorMsg("Error in SusyResonanceWidths::allowcalc: "
    "unable to reset decay table.", idStream.str(), true);
  return done;
}

//==========================================================================

// The ResonanceSquark class
// Derived class for Squark resonances

//--------------------------------------------------------------------------

// Set up channels

bool ResonanceSquark::getChannels(int idPDG){

  idPDG = abs(idPDG);

  int ksusy = 1000000;
  if (idPDG < ksusy) return false;
  if(idPDG % ksusy >= 7 || idPDG % ksusy < 1) return false;

  ParticleDataEntry* squarkEntryPtr
    = particleDataPtr->particleDataEntryPtr(idPDG);

  // Delete any decay channels read
  squarkEntryPtr->clearChannels();

  if(idPDG % 2 == 0) { // up-type squarks

    /*
    // Gravitino
    squarkEntryPtr->addChannel(1, 0.0, 103, 1000039, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000039, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000039, 6);
    */

    // Gaugino - quark
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000024, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000024, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000037, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000037, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000037, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000022, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000022, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000022, 6);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000023, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000023, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000023, 6);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000025, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000025, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000025, 6);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000035, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000035, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000035, 6);

    // Squark - Gauge boson
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000001, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000003, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000005, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000001, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000003, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000005, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000001, -37);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000003, -37);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000005, -37);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000001, -37);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000003, -37);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000005, -37);

    // Gluino - quark
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000021, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000021, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000021, 6);

    // lepton - quark via LQD
    squarkEntryPtr->addChannel(1, 0.0, 0, -11, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, -11, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, -11, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, -13, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, -13, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, -13, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, -15, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, -15, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, -15, 5);

    // quark - quark via UDD
    squarkEntryPtr->addChannel(1, 0.0, 0, -1 ,-3);
    squarkEntryPtr->addChannel(1, 0.0, 0, -1 ,-5);
    squarkEntryPtr->addChannel(1, 0.0, 0, -3 ,-5);


  } else { //down-type squarks

    // Gaugino - quark
    squarkEntryPtr->addChannel(1, 0.0, 0, -1000024, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, -1000037, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, -1000024, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, -1000037, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, -1000024, 6);
    squarkEntryPtr->addChannel(1, 0.0, 0, -1000037, 6);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000022, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000022, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000022,  5);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000023, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000023, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000023, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000025, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000025, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000025, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000035, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000035, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000035, 5);

    // Squark - Gauge boson
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000002, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000004, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000006, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000002, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000004, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000006, -24);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000002, -37);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000004, -37);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000006, -37);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000002, -37);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000004, -37);
    squarkEntryPtr->addChannel(1, 0.0, 0, 2000006, -37);

    // Gluino - quark
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000021, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000021, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, 1000021, 5);

    // lepton - quark via LQD
    squarkEntryPtr->addChannel(1, 0.0, 0, -12, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, -12, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, -12, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, -14, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, -14, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, -14, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, -16, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, -16, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, -16, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, 12 ,1);
    squarkEntryPtr->addChannel(1, 0.0, 0, 11 ,2);
    squarkEntryPtr->addChannel(1, 0.0, 0, 12, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, 11, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, 12, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, 11, 6);
    squarkEntryPtr->addChannel(1, 0.0, 0, 14, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, 13, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, 14, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, 13, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, 14, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, 13, 6);
    squarkEntryPtr->addChannel(1, 0.0, 0, 16, 1);
    squarkEntryPtr->addChannel(1, 0.0, 0, 15, 2);
    squarkEntryPtr->addChannel(1, 0.0, 0, 16, 3);
    squarkEntryPtr->addChannel(1, 0.0, 0, 15, 4);
    squarkEntryPtr->addChannel(1, 0.0, 0, 16, 5);
    squarkEntryPtr->addChannel(1, 0.0, 0, 15, 6);


    // quark - quark via UDD
    squarkEntryPtr->addChannel(1, 0.0, 0, -2, -1);
    squarkEntryPtr->addChannel(1, 0.0, 0, -2, -3);
    squarkEntryPtr->addChannel(1, 0.0, 0, -2, -5);
    squarkEntryPtr->addChannel(1, 0.0, 0, -4, -1);
    squarkEntryPtr->addChannel(1, 0.0, 0, -4, -3);
    squarkEntryPtr->addChannel(1, 0.0, 0, -4, -5);
    squarkEntryPtr->addChannel(1, 0.0, 0, -6, -1);
    squarkEntryPtr->addChannel(1, 0.0, 0, -6, -3);
    squarkEntryPtr->addChannel(1, 0.0, 0, -6, -5);
  }

  return true;

}

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceSquark::initConstants() {

  // Locally stored properties and couplings.
  s2W = coupSUSYPtr->sin2W;

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceSquark::calcPreFac(bool) {

  // Common coupling factors.
  alpS   = coupSUSYPtr->alphaS(mHat * mHat );
  alpEM  = coupSUSYPtr->alphaEM(mHat * mHat);
  preFac = 1.0 / (s2W * pow(mHat,3));
  ps *= mHat * mHat;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceSquark::calcWidth(bool) {

  // Squark type -- in u_i/d_i and generation
  int ksusy = 1000000;
  bool idown = (abs(idRes)%2 == 0 ? false : true);
  int isq = (abs(idRes)/ksusy == 2) ?
    (abs(idRes)%10+1)/2 + 3: (abs(idRes)%10+1)/2;
  // int isqgen = (abs(idRes)%10 + 1)/2;

  // Check that mass is above threshold.
  if (ps == 0.) return;
  else{
    // Two-body decays
    kinFac = (mHat * mHat - mf1 * mf1 - mf2 * mf2);

    double fac = 0.0 , wid = 0.0;

    //RPV decays
    //Case 1a:  UDD-type
    if (id1Abs < 7 && id2Abs < 7){

      // Quark generations
      int iq1 = (id1Abs + 1)/2;
      int iq2 = (id2Abs + 1)/2;

      // Check for RPV UDD couplings
      if (!coupSUSYPtr->isUDD) {widNow = 0; return;}

      // ~q -> q_i + q_j

      fac = 2.0 * kinFac / (16.0 * M_PI * pow(mHat,3));
      wid = 0.0;
      if (idown) {
        if ((id1Abs+id2Abs)%2 == 1){
          if (id1Abs%2==1)
            for (int isq2 = 1; isq2 < 4; isq2++)
              wid += norm(coupSUSYPtr->rvUDD[iq2][iq1][isq2]
                   * coupSUSYPtr->Rdsq[isq][isq2+3]);
          else
            for (int isq2 = 1; isq2 < 4; isq2++)
              wid += norm(coupSUSYPtr->rvUDD[iq1][iq2][isq2]
                   * coupSUSYPtr->Rdsq[isq][isq2+3]);
        }
      }
      else {
        if ((id1Abs+id2Abs)%2 != 0) widNow = 0.0;
        else
          for (int isq2 = 1; isq2 < 4; isq2++)
            wid += norm(coupSUSYPtr->rvUDD[isq2][iq1][iq2]
                 * coupSUSYPtr->Rusq[isq][isq2+3]);
      }
  }

    //Case 1b:  LQD-type
    else if (id1Abs < 17 && id2Abs < 7){
      if (!coupSUSYPtr->isLQD) {widNow = 0; return;}

      int ilep = (id1Abs - 9)/2;
      int iq = (id2Abs + 1)/2;

      fac = kinFac / (16.0 * M_PI * pow(mHat,3));
      wid = 0.0;
      if (idown){
        if (iq%2 == 0){
          // q is up-type; ~q is right-handed down type
          for (int isq2=1; isq2<3; isq2++)
            wid += norm(coupSUSYPtr->Rdsq[isq][isq2+3]
                 * coupSUSYPtr->rvLQD[ilep][iq][isq2]);
        }else{
          //q is down type; ~q left-handed down-type
          for (int isq2=1; isq2<3; isq2++)
            wid += norm(coupSUSYPtr->Rdsq[isq][isq2]
                 * coupSUSYPtr->rvLQD[ilep][isq2][isq2]);
        }
      }
      else{
        if (iq%2 == 0) {widNow = 0.0; return;}
        // q is down type; ~q is left-handed up-type
        for (int isq2=1; isq2<3; isq2++)
          wid += norm(coupSUSYPtr->Rusq[isq][isq2]
               * coupSUSYPtr->rvLQD[ilep][isq2][iq]);
      }
    }

    //Case 2: quark + gaugino
    else if (id1Abs > ksusy && id2Abs < 7) {

      int iq = (id2Abs + 1)/2;

      // ~q -> ~g + q
      if (id1Abs == 1000021 && idRes%10 == id2Abs) {
        // Removed factor of s2W in denominator: strong process -- no EW
        fac = 2.0 * alpS / (3.0 * pow3(mHat));
        if (idown)
          wid = kinFac * (norm(coupSUSYPtr->LsddG[isq][iq])
              + norm(coupSUSYPtr->RsddG[isq][iq]))
              - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsddG[isq][iq]
              * conj(coupSUSYPtr->RsddG[isq][iq]));
        else
          wid = kinFac * (norm(coupSUSYPtr->LsuuG[isq][iq])
              + norm(coupSUSYPtr->RsuuG[isq][iq]))
              - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsuuG[isq][iq]
              * conj(coupSUSYPtr->RsuuG[isq][iq]));
      }
      else
        for (int i=1; i<6 ; i++){
          // ~q -> ~chi0 + q
          if (coupSUSYPtr->idNeut(i)==id1Abs && idRes%2 == id2Abs%2){
            fac = alpEM *  preFac / (2.0 * (1 - s2W));
            if (idown)
              wid = kinFac * (norm(coupSUSYPtr->LsddX[isq][iq][i])
                  + norm(coupSUSYPtr->RsddX[isq][iq][i]))
                  - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsddX[isq][iq][i]
                  * conj(coupSUSYPtr->RsddX[isq][iq][i]));
            else
              wid = kinFac * (norm(coupSUSYPtr->LsuuX[isq][iq][i])
                  + norm(coupSUSYPtr->RsuuX[isq][iq][i]))
                  - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsuuX[isq][iq][i]
                  * conj(coupSUSYPtr->RsuuX[isq][iq][i]));
          }

          // ~q -> chi- + q
          else if (i < 3 && coupSUSYPtr->idChar(i)==id1Abs
            && idRes%2 != id2Abs%2){

            fac = alpEM *  preFac / (4.0 * (1 - s2W));
            if (idown)
              wid = kinFac * (norm(coupSUSYPtr->LsduX[isq][iq][i])
                  + norm(coupSUSYPtr->RsduX[isq][iq][i]))
                  - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsduX[isq][iq][i]
                  * conj(coupSUSYPtr->RsduX[isq][iq][i]));
            else
              wid = kinFac * (norm(coupSUSYPtr->LsudX[isq][iq][i])
                  + norm(coupSUSYPtr->RsudX[isq][iq][i]))
                  - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsudX[isq][iq][i]
                  * conj(coupSUSYPtr->RsudX[isq][iq][i]));
          }
        }
    }

    //Case 3: ~q_i -> ~q_j + Z/W
    else if (id1Abs > ksusy && id1Abs%100 < 7
      && (id2Abs == 23 || id2Abs == 24)){

      // factor of lambda^(3/2) = ps^(3/2) ;
      fac = alpEM * preFac/(16.0 * pow2(particleDataPtr->m0(id2Abs))
          * (1.0 - s2W)) * pow2(ps) ;

      int isq2 = (id1Abs/ksusy == 2) ? (id1Abs%10+1)/2 + 3: (id1Abs%10+1)/2;

      if (id2Abs == 23 && id1Abs%2 == idRes%2){
        if (idown)
          wid = norm(coupSUSYPtr->LsdsdZ[isq][isq2]
                     + coupSUSYPtr->RsdsdZ[isq][isq2]);
        else
          wid = norm(coupSUSYPtr->LsusuZ[isq][isq2]
                     + coupSUSYPtr->RsusuZ[isq][isq2]);
      }
      else if (id2Abs == 24 && id1Abs%2 != idRes%2){
        if (idown)
          wid = norm(coupSUSYPtr->LsusdW[isq2][isq]);
        else
          wid = norm(coupSUSYPtr->LsusdW[isq][isq2]);
      }
    }

    // TODO: Case ~q_i -> ~q_j + h/H
    widNow = fac * wid * ps * pow2(mHat);
    if (DBSUSY) cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs
                  <<" Width: "<<widNow<<endl;
    return;
  }

}

//==========================================================================

// The ResonanceGluino class
// Derived class for Gluino resonances

//--------------------------------------------------------------------------

// Set up Channels

bool ResonanceGluino::getChannels(int idPDG){

  idPDG = abs(idPDG);
  if (idPDG != 1000021) return false;

  ParticleDataEntry* gluinoEntryPtr
    = particleDataPtr->particleDataEntryPtr(idPDG);

  // Delete any decay channels read
  gluinoEntryPtr->clearChannels();

  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000001, -1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000001, 1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000001 ,-3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000001, 3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000001 ,-5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000001, 5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000001 ,-1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000001, 1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000001 ,-3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000001, 3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000001 ,-5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000001, 5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000002 ,-2);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000002, 2);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000002 ,-4);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000002, 4);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000002 ,-6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000002, 6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000002 ,-2);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000002, 2);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000002 ,-4);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000002, 4);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000002 ,-6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000002, 6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000003 ,-1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000003, 1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000003 ,-3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000003, 3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000003 ,-5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000003, 5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000003 ,-1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000003, 1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000003 ,-3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000003, 3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000003 ,-5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000003, 5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000004 ,-2);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000004, 2);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000004 ,-4);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000004, 4);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000004 ,-6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000004, 6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000004 ,-2);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000004, 2);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000004 ,-4);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000004, 4);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000004 ,-6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000004, 6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000005 ,-1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000005, 1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000005 ,-3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000005, 3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000005 ,-5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000005, 5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000005 ,-1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000005, 1);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000005 ,-3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000005, 3);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000005 ,-5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000005, 5);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000006 ,-6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000006, 6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000006 ,-2);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000006, 2);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 1000006 ,-4);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -1000006, 4);
  gluinoEntryPtr->addChannel(1, 0.0, 0, 2000006 ,-6);
  gluinoEntryPtr->addChannel(1, 0.0, 0, -2000006, 6);

  return true;
}

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceGluino::initConstants() {
  return;
}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceGluino::calcPreFac(bool) {

  // Common coupling factors.
  alpS  = coupSUSYPtr->alphaS(mHat * mHat);
  preFac = alpS /( 8.0 * pow(mHat,3));
  return;

}


//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceGluino::calcWidth(bool) {

  widNow = 0.0;
  if (ps == 0.) return;
  kinFac = (mHat * mHat - mf1 * mf1 + mf2 * mf2);

  if (id1Abs > 1000000 && (id1Abs % 100) < 7 && id2Abs < 7) {

    int isq = (abs(id1Abs)/1000000 == 2) ? (abs(id1Abs)%10+1)/2 + 3
                                         : (abs(id1Abs)%10+1)/2;
    bool idown = id2Abs%2;
    int iq = (id2Abs + 1)/2;

    // ~g -> ~q + q
    if (idown){
      widNow = kinFac * (norm(coupSUSYPtr->LsddG[isq][iq])
             + norm(coupSUSYPtr->RsddG[isq][iq]))
             + 4.0 * mHat * mf2 * real(coupSUSYPtr->LsddG[isq][iq]
             * conj(coupSUSYPtr->RsddG[isq][iq]));
    }
    else{
      widNow = kinFac * (norm(coupSUSYPtr->LsuuG[isq][iq])
             + norm(coupSUSYPtr->RsuuG[isq][iq]))
             + 4.0 * mHat * mf2 * real(coupSUSYPtr->LsuuG[isq][iq]
             * conj(coupSUSYPtr->RsuuG[isq][iq]));
    }
    widNow = widNow * preFac * ps * pow2(mHat);
    if (DBSUSY) {
      cout<<"Gluino:: id1:"<<id1Abs<<" id2:"<<id2Abs<<" Width: ";
      cout<<scientific<<widNow<<endl;
    }
    return;
  }
}

//==========================================================================

//  Class ResonanceNeut
//  Derived class for Neutralino Resonances
//
//--------------------------------------------------------------------------

// Set up Channels

bool ResonanceNeut::getChannels(int idPDG){

  idPDG = abs(idPDG);

  int iNeut = coupSUSYPtr->typeNeut(idPDG);
  if (iNeut < 1) return false;

  ParticleDataEntry* neutEntryPtr
    = particleDataPtr->particleDataEntryPtr(idPDG);

  // Delete any decay channels read
  neutEntryPtr->clearChannels();

  // RPV channels

  neutEntryPtr->addChannel(1, 0.0, 0, -12, -13, 11);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 13, -11);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -13, 13);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 13, -13);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -13, 15);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 13, -15);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -15, 11);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 15, -11);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -15, 13);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 15, -13);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -15, 15);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 15, -15);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -11, 11);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 11, -11);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -11, 13);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 11, -13);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -11, 15);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 11, -15);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -15, 11);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 15, -11);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -15, 13);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 15, -13);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -15, 15);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 15, -15);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -11, 11);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 11, -11);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -11, 13);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 11, -13);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -11, 15);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 11, -15);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -13, 11);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 13, -11);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -13, 13);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 13, -13);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -13, 15);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 13, -15);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -1, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 1, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -11, -2, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 11, 2, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -1, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 1, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -11, -2, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 11, 2, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -1, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 1, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -11, -2, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 11, 2, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -3, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 3, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -11, -4, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 11, 4, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -3, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 3, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -11, -4, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 11, 4, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -3, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 3, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -11, -4, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 11, 4, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -5, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 5, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -11, -6, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 11, 6, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -12, -5, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 5, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -11, -6, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 11, 6, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, 12 ,-5 ,5);
  neutEntryPtr->addChannel(1, 0.0, 0, 12, 5, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -11, -6, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 11, 6, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -1, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 1, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -13, -2, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 13, 2, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -1, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 1, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -13, -2, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 13, 2, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -1, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 1, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -13, -2, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 13, 2, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -3, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 3, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -13, -4, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 13, 4, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -3, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 3, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -13, -4, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 13, 4, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -3, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 3, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -13, -4, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 13, 4, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -5, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 5, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -13, -6, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 13, 6, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -5, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 5, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -13, -6, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 13, 6, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -14, -5, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 14, 5, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -13, -6, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 13, 6, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -1, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 1, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -15, -2, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 15, 2, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -1, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 1, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -15, -2, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 15, 2, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -1, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 1, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -15, -2, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 15, 2, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -3, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 3, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -15, -4, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 15, 4, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -3, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 3, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -15, -4, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 15, 4, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -3, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 3, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -15, -4, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 15, 4, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -5, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 5, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -15, -6, 1);
  neutEntryPtr->addChannel(1, 0.0, 0, 15, 6, -1);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -5, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 5, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -15, -6, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, 15, 6, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, -16, -5, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 16, 5, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -15, -6, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, 15, 6, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, -2 ,-1 ,-3);
  neutEntryPtr->addChannel(1, 0.0, 0, 2, 1, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, -2, -1, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, 2, 1, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, -2, -3, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, 2, 3, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, -4, -1, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, 4, 1, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, -4, -1, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, 4, 1, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, -4, -3, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, 4, 3, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, -6, -1, -3);
  neutEntryPtr->addChannel(1, 0.0, 0, 6, 1, 3);
  neutEntryPtr->addChannel(1, 0.0, 0, -6, -1, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, 6, 1, 5);
  neutEntryPtr->addChannel(1, 0.0, 0, -6, -3, -5);
  neutEntryPtr->addChannel(1, 0.0, 0, 6, 3, 5);

  if (iNeut > 1) {

    neutEntryPtr->addChannel(1, 0.0, 0, 1000022, 22);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000022, 23);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000022, 25);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000022, 35);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000022, 36);

    if ( iNeut > 2) {
    neutEntryPtr->addChannel(1, 0.0, 0, 1000023, 22);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000023, 23);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000023, 25);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000023, 35);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000023, 36);
    }

    if (iNeut > 3) {
    neutEntryPtr->addChannel(1, 0.0, 0, 1000025, 22);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000025, 23);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000025, 25);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000025, 35);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000025, 36);
    }

    if (iNeut > 4) {
    neutEntryPtr->addChannel(1, 0.0, 0, 1000035, 22);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000035, 23);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000035, 25);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000035, 35);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000035, 36);
    }

    neutEntryPtr->addChannel(1, 0.0, 0, 1000024, -24);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000024, 24);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000037, -24);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000037, 24);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000024, -37);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000024, 37);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000037, -37);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000037, 37);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000011, -11);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000011, 11);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000011, -11);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000011, 11);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000012, -12);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000012, 12);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000013, -13);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000013, 13);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000013, -13);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000013, 13);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000014, -14);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000014, 14);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000015, -15);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000015, 15);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000015, -15);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000015, 15);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000016, -16);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000016, 16);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000001, -1);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000001, 1);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000001, -3);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000001, 3);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000001, -5);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000001, 5);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000001, -1);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000001, 1);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000001, -3);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000001, 3);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000001, -5);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000001, 5);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000002, -2);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000002, 2);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000002, -4);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000002, 4);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000002, -6);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000002, 6);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000002, -2);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000002, 2);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000002, -4);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000002, 4);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000002, -6);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000002, 6);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000003, -1);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000003, 1);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000003, -3);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000003, 3);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000003, -5);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000003, 5);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000003, -1);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000003, 1);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000003, -3);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000003, 3);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000003, -5);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000003, 5);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000004, -2);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000004, 2);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000004, -4);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000004, 4);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000004, -6);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000004, 6);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000004, -2);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000004, 2);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000004, -4);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000004, 4);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000004, -6);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000004, 6);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000005, -1);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000005, 1);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000005, -3);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000005, 3);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000005, -5);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000005, 5);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000005, -1);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000005, 1);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000005, -3);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000005, 3);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000005, -5);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000005, 5);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000006, -6);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000006, 6);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000006, -2);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000006, 2);
    neutEntryPtr->addChannel(1, 0.0, 0, 1000006, -4);
    neutEntryPtr->addChannel(1, 0.0, 0, -1000006, 4);
    neutEntryPtr->addChannel(1, 0.0, 0, 2000006, -6);
    neutEntryPtr->addChannel(1, 0.0, 0, -2000006, 6);

    // Modes involving right-handed sneutrinos are not included by default,
    // but can be added by hand, by uncommenting the following lines.
    // neutEntryPtr->addChannel(1, 0.0, 0, 2000012, -12);
    // neutEntryPtr->addChannel(1, 0.0, 0, -2000012, 12);
    // neutEntryPtr->addChannel(1, 0.0, 0, 2000014, -14);
    // neutEntryPtr->addChannel(1, 0.0, 0, -2000014, 14);
    // neutEntryPtr->addChannel(1, 0.0, 0, 2000016, -16);
    // neutEntryPtr->addChannel(1, 0.0, 0, -2000016, 16);



  }
  return true;
}

//--------------------------------------------------------------------------

void ResonanceNeut::initConstants() {

  s2W = coupSUSYPtr->sin2W;

  // Initialize functions for calculating 3-body widths
  // psi.setPointers(particleDataPtr,coupSUSYPtr,infoPtr);
  // phi.setPointers(particleDataPtr,coupSUSYPtr,infoPtr);
  // upsil.setPointers(particleDataPtr,coupSUSYPtr,infoPtr);

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void  ResonanceNeut::calcPreFac(bool) {

  // Common coupling factors.
  alpEM  = coupSUSYPtr->alphaEM(mHat * mHat);
  preFac = alpEM / (8.0 * s2W * pow(mHat,3));
  return;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.
void  ResonanceNeut::calcWidth(bool){

  widNow = 0.0;

  if (ps == 0.) return;
  else if (mult ==2){
    // Two-body decays

    kinFac = mHat * mHat - mf1 * mf1 + mf2 * mf2;
    kinFac2 = pow(mHat,4) + pow(mf1,4) - 2.0 * pow(mf2,4)
            + pow2(mHat) * pow2(mf2) + pow2(mf1) * pow2(mf2)
            - 2.0 * pow2(mHat) * pow2(mf1);

    // Stable lightest neutralino
    if (idRes == 1000022) return;

    double fac = 0.0;
    int iNeut1 = coupSUSYPtr->typeNeut(idRes);
    int iNeut2 = coupSUSYPtr->typeNeut(id1Abs);
    int iChar1 = coupSUSYPtr->typeChar(id1Abs);

    if (iNeut2>0 && id2Abs == 23){
      // ~chi0_i -> chi0_j + Z
      fac = kinFac2 * (norm(coupSUSYPtr->OLpp[iNeut1][iNeut2])
          + norm(coupSUSYPtr->ORpp[iNeut1][iNeut2]));
      fac -= 12.0 * mHat * mf1 * pow2(mf2)
           * real(coupSUSYPtr->OLpp[iNeut1][iNeut2]
           * conj(coupSUSYPtr->ORpp[iNeut1][iNeut2]));
      fac /= pow2(mf2) * (1.0 - s2W);
    }
    else if (iChar1>0 && id2Abs==24){
      // ~chi0_i -> chi+_j + W- (or c.c.)

      fac = kinFac2 * (norm(coupSUSYPtr->OL[iNeut1][iChar1])
          + norm(coupSUSYPtr->OR[iNeut1][iChar1]));
      fac -= 12.0 * mHat * mf1 * pow2(mf2)
           * real(coupSUSYPtr->OL[iNeut1][iChar1]
           * conj(coupSUSYPtr->OR[iNeut1][iChar1]));
      fac /= pow2(mf2);
    }
    else if (id1Abs > 1000000 && id1Abs%100 < 7 && id2Abs < 7){
      // ~chi0_k -> ~q + q
      bool idown = (id1Abs%2 == 1);
      int iq = (id2Abs + 1 )/ 2;
      int isq = (abs(id1Abs)/1000000 == 2) ? (abs(id1Abs)%10+1)/2 + 3
                                           : (abs(id1Abs)%10+1)/2;

      if (idown){
        fac  = kinFac * (norm(coupSUSYPtr->LsddX[isq][iq][iNeut1])
             + norm(coupSUSYPtr->RsddX[isq][iq][iNeut1]));
        fac += 4.0 * mHat * mf2 * real(coupSUSYPtr->LsddX[isq][iq][iNeut1]
             * conj(coupSUSYPtr->RsddX[isq][iq][iNeut1]));
      }
      else{
        fac = kinFac * (norm(coupSUSYPtr->LsuuX[isq][iq][iNeut1])
            + norm(coupSUSYPtr->RsuuX[isq][iq][iNeut1]));
        fac += 4.0 * mHat * mf2 * real(coupSUSYPtr->LsuuX[isq][iq][iNeut1]
             * conj(coupSUSYPtr->RsuuX[isq][iq][iNeut1]));
      }
      // Extra multiplicative factor of 3 over sleptons
      fac *= 6.0/(1 - s2W);
    }
    else if (id1Abs > 2000010 && id1Abs%2 == 0 ) {
      // Check for right-handed neutralinos.
      widNow = 0;
    }
    else if (id1Abs > 1000000 && id1Abs%100 > 10 && id1Abs%100 < 17
      && id2Abs < 17){

      // ~chi0_k -> ~l + l
      bool idown = id2Abs%2;
      int il = (id2Abs - 9)/ 2;
      int isl = (abs(id1Abs)/1000000 == 2) ? (abs(id1Abs)%10+1)/2 + 3
                                           : (abs(id1Abs)%10+1)/2;

      if (idown){
        fac  = kinFac * (norm(coupSUSYPtr->LsllX[isl][il][iNeut1])
             + norm(coupSUSYPtr->RsllX[isl][il][iNeut1]));
        fac += 4.0 * mHat * mf2 * real(coupSUSYPtr->LsllX[isl][il][iNeut1]
             * conj(coupSUSYPtr->RsllX[isl][il][iNeut1]));
      }
      else{
        fac = kinFac * (norm(coupSUSYPtr->LsvvX[isl][il][iNeut1]));
      }
      fac *= 2.0/(1 - s2W);
    }
    // TODO: Decays in higgs
    // Final width for 2-body decays
    widNow = fac * preFac * ps * pow2(mHat);
    if (DBSUSY) {
      cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs<<" Width: ";
      cout<<scientific<<widNow<<endl;
    }
  }
  else {
    //RPV 3-body decays

    //Case: UDD-type (TO BE RE-DONE to fix bug)
    return;

    }

  // Normalisation.  Extra factor of 2 to account for the fact that
  // d_i, d_j will always be ordered in ascending order. N_c! = 6
  widNow *= 12.0 /(pow3(2.0 * M_PI * mHat) * 32.0);

  return;
}

//==========================================================================

//  Class ResonanceChar
//  Derived class for Neutralino Resonances
//  Decays into higgses/sleptons not yet implemented

//--------------------------------------------------------------------------

// Set up Channels

bool ResonanceChar::getChannels(int idPDG){

  idPDG = abs(idPDG);
  int iChar = coupSUSYPtr->typeChar(idPDG);
  if (iChar < 1) return false;

  ParticleDataEntry* charEntryPtr
    = particleDataPtr->particleDataEntryPtr(idPDG);

  // Delete any decay channels read
  charEntryPtr->clearChannels();

  charEntryPtr->addChannel(1, 0.0, 0, 1000022, 24);
  charEntryPtr->addChannel(1, 0.0, 0, 1000023, 24);
  charEntryPtr->addChannel(1, 0.0, 0, 1000025, 24);
  charEntryPtr->addChannel(1, 0.0, 0, 1000035, 24);
  charEntryPtr->addChannel(1, 0.0, 0, 1000022, 37);
  charEntryPtr->addChannel(1, 0.0, 0, 1000023, 37);
  charEntryPtr->addChannel(1, 0.0, 0, 1000025, 37);
  charEntryPtr->addChannel(1, 0.0, 0, 1000035, 37);
  charEntryPtr->addChannel(1, 0.0, 0, 1000012, -11);
  charEntryPtr->addChannel(1, 0.0, 0, -1000011, 12);
  charEntryPtr->addChannel(1, 0.0, 0, -2000011, 12);
  charEntryPtr->addChannel(1, 0.0, 0, 1000014, -13);
  charEntryPtr->addChannel(1, 0.0, 0, -1000013, 14);
  charEntryPtr->addChannel(1, 0.0, 0, -2000013, 14);
  charEntryPtr->addChannel(1, 0.0, 0, 1000016, -15);
  charEntryPtr->addChannel(1, 0.0, 0, -1000015, 16);
  charEntryPtr->addChannel(1, 0.0, 0, -2000015, 16);
  charEntryPtr->addChannel(1, 0.0, 0, 1000002, -1);
  charEntryPtr->addChannel(1, 0.0, 0, 1000002, -3);
  charEntryPtr->addChannel(1, 0.0, 0, 1000002, -5);
  charEntryPtr->addChannel(1, 0.0, 0, 2000002, -1);
  charEntryPtr->addChannel(1, 0.0, 0, 2000002, -3);
  charEntryPtr->addChannel(1, 0.0, 0, 2000002, -5);
  charEntryPtr->addChannel(1, 0.0, 0, -1000001, 2);
  charEntryPtr->addChannel(1, 0.0, 0, -1000001, 4);
  charEntryPtr->addChannel(1, 0.0, 0, -1000001, 6);
  charEntryPtr->addChannel(1, 0.0, 0, -2000001, 2);
  charEntryPtr->addChannel(1, 0.0, 0, -2000001, 4);
  charEntryPtr->addChannel(1, 0.0, 0, -2000001, 6);
  charEntryPtr->addChannel(1, 0.0, 0, 1000004, -1);
  charEntryPtr->addChannel(1, 0.0, 0, 1000004, -3);
  charEntryPtr->addChannel(1, 0.0, 0, 1000004, -5);
  charEntryPtr->addChannel(1, 0.0, 0, 2000004, -1);
  charEntryPtr->addChannel(1, 0.0, 0, 2000004, -3);
  charEntryPtr->addChannel(1, 0.0, 0, 2000004, -5);
  charEntryPtr->addChannel(1, 0.0, 0, -1000003, 2);
  charEntryPtr->addChannel(1, 0.0, 0, -1000003, 4);
  charEntryPtr->addChannel(1, 0.0, 0, -1000003, 6);
  charEntryPtr->addChannel(1, 0.0, 0, -2000003, 2);
  charEntryPtr->addChannel(1, 0.0, 0, -2000003, 4);
  charEntryPtr->addChannel(1, 0.0, 0, -2000003, 6);
  charEntryPtr->addChannel(1, 0.0, 0, 1000006, -1);
  charEntryPtr->addChannel(1, 0.0, 0, 1000006, -3);
  charEntryPtr->addChannel(1, 0.0, 0, 1000006, -5);
  charEntryPtr->addChannel(1, 0.0, 0, 2000006, -1);
  charEntryPtr->addChannel(1, 0.0, 0, 2000006, -3);
  charEntryPtr->addChannel(1, 0.0, 0, 2000006, -5);
  charEntryPtr->addChannel(1, 0.0, 0, -1000005, 2);
  charEntryPtr->addChannel(1, 0.0, 0, -1000005, 4);
  charEntryPtr->addChannel(1, 0.0, 0, -1000005, 6);
  charEntryPtr->addChannel(1, 0.0, 0, -2000005, 2);
  charEntryPtr->addChannel(1, 0.0, 0, -2000005, 4);
  charEntryPtr->addChannel(1, 0.0, 0, -2000005, 6);

  if (iChar > 1) {
    charEntryPtr->addChannel(1, 0.0, 0, 1000024, 23);
    charEntryPtr->addChannel(1, 0.0, 0, 1000024, 25);
    charEntryPtr->addChannel(1, 0.0, 0, 1000024, 35);
    charEntryPtr->addChannel(1, 0.0, 0, 1000024, 36);
  }

  // Modes involving right-handed sneutrinos are not included by default,
  // but can be added by hand, by uncommenting the following lines.
  // charEntryPtr->addChannel(1, 0.0, 0, 2000012, -11);
  // charEntryPtr->addChannel(1, 0.0, 0, 2000014, -13);
  // charEntryPtr->addChannel(1, 0.0, 0, 2000016, -15);

  return true;

}

//--------------------------------------------------------------------------

void ResonanceChar::initConstants() {

  s2W = coupSUSYPtr->sin2W;
  return;

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void  ResonanceChar::calcPreFac(bool) {

  alpEM  = coupSUSYPtr->alphaEM(mHat * mHat);
  preFac = alpEM / (8.0 * s2W * pow(mHat,3));
  return;

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void  ResonanceChar::calcWidth(bool) {

  widNow = 0.0;
  if (ps == 0.) return;

  if (mult ==2){
    double fac = 0.0;
    kinFac = mHat * mHat - mf1 * mf1 + mf2 * mf2;
    kinFac2 = pow(mHat,4) + pow(mf1,4) - 2.0 * pow(mf2,4)
      + pow2(mHat) * pow2(mf2) + pow2(mf1)
      * pow2(mf2) - 2.0 * pow2(mHat) * pow2(mf1);

    int idChar1 = coupSUSYPtr->typeChar(idRes);
    int idChar2 = coupSUSYPtr->typeChar(id1Abs);
    int idNeut1 = coupSUSYPtr->typeNeut(id1Abs);

    if (idChar2>0 && id2Abs == 23){
      // ~chi_i -> chi_j + Z
      fac = kinFac2 * (norm(coupSUSYPtr->OLp[idChar1][idChar2])
          + norm(coupSUSYPtr->ORp[idChar1][idChar2]));
      fac -= 12.0 * mHat * mf1 * pow2(mf2)
           * real(coupSUSYPtr->OLp[idChar1][idChar2]
           * conj(coupSUSYPtr->ORp[idChar1][idChar2]));
      fac /= pow2(mf2) * (1.0 - s2W);
    }
    else if (idNeut1>0 && id2Abs==24){
      // ~chi_i -> chi0_j + W- (or c.c.)

      fac  = kinFac2 * (norm(coupSUSYPtr->OL[idNeut1][idChar1])
           + norm(coupSUSYPtr->OR[idNeut1][idChar1]));
      fac -= 12.0 * mHat * mf1 * pow2(mf2)
           * real(coupSUSYPtr->OL[idNeut1][idChar1]
           * conj(coupSUSYPtr->OR[idNeut1][idChar1]));
      fac /= pow2(mf2);
    }
    else if (id1Abs > 1000000 && id1Abs%100 < 7 && id2Abs < 7){
      // ~chi0_k -> ~q + q
      bool idown = (id1Abs%2 == 1);
      int iq = (id2Abs + 1 )/ 2;
      int isq = (abs(id1Abs)/1000000 == 2) ? (abs(id1Abs)%10+1)/2 + 3
                                           : (abs(id1Abs)%10+1)/2;

      if (idown){
        fac  = kinFac * (norm(coupSUSYPtr->LsduX[isq][iq][idChar1])
             + norm(coupSUSYPtr->RsduX[isq][iq][idChar1]));
        fac += 4.0 * mHat * mf2
             * real(coupSUSYPtr->LsduX[isq][iq][idChar1]
             * conj(coupSUSYPtr->RsduX[isq][iq][idChar1]));
      }
      else{
        fac  = kinFac * (norm(coupSUSYPtr->LsudX[isq][iq][idChar1])
             + norm(coupSUSYPtr->RsudX[isq][iq][idChar1]));
        fac += 4.0 * mHat * mf2
             * real(coupSUSYPtr->LsudX[isq][iq][idChar1]
             * conj(coupSUSYPtr->RsudX[isq][iq][idChar1]));
      }
      fac *= 6.0/(1 - s2W);
    }
    else if (id1Abs > 2000010 && id1Abs%2 == 0 ) {
      // Check for right-handed neutralinos.
      widNow = 0;
    }
    else if (id1Abs > 1000000 && id1Abs%100 > 10 && id1Abs%100 < 17
      && id2Abs < 17){
      // ~chi+_k -> ~l + l
      bool idown = id2Abs%2;
      int il = (id2Abs - 9)/ 2;
      int isl = (abs(id1Abs)/1000000 == 2) ? (abs(id1Abs)%10+1)/2 + 3
                                           : (abs(id1Abs)%10+1)/2;

      if (idown){
        fac  = kinFac * (norm(coupSUSYPtr->LslvX[isl][il][idChar1])
             + norm(coupSUSYPtr->RslvX[isl][il][idChar1]));
        fac += 4.0 * mHat * mf2 * real(coupSUSYPtr->LslvX[isl][il][idChar1]
             * conj(coupSUSYPtr->RslvX[isl][il][idChar1]));
      }
      else{
        fac = kinFac * (norm(coupSUSYPtr->LsvlX[isl][il][idChar1]));
      }
      fac *= 2.0/(1 - s2W);
    }

    // TODO: Decays in higgs
    widNow = fac * preFac * ps * pow2(mHat);
    if (DBSUSY) {
      cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs<<" Width: ";
      cout<<scientific<<widNow<<endl;
    }
  }else{
    //TODO: Implement Chargino 3-body decays
  }
  return;
}

//==========================================================================
// The ResonanceSlepton class
// Derived class for Slepton (and sneutrino) resonances

//--------------------------------------------------------------------------

// Set up Channels

bool ResonanceSlepton::getChannels(int idPDG){

  idPDG = abs(idPDG);

  int ksusy = 1000000;
  if (idPDG < ksusy) return false;
  if(idPDG % ksusy < 7 || idPDG % ksusy > 17) return false;

  ParticleDataEntry* slepEntryPtr
    = particleDataPtr->particleDataEntryPtr(idPDG);

  // Delete any decay channels read
  slepEntryPtr->clearChannels();

  if(idPDG % 2 == 1) {

    slepEntryPtr->addChannel(1, 0.0, 0, -1000024, 16);
    slepEntryPtr->addChannel(1, 0.0, 0, -1000037, 16);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000022, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000023, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000025, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000035, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000016, -24);
    slepEntryPtr->addChannel(1, 0.0, 0, 2000016, -24);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000016, -37);
    slepEntryPtr->addChannel(1, 0.0, 0, 2000016, -37);
    slepEntryPtr->addChannel(1, 0.0, 0, 12, 13);
    slepEntryPtr->addChannel(1, 0.0, 0, 12, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, 14, 11);
    slepEntryPtr->addChannel(1, 0.0, 0, 14, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, 16, 11);
    slepEntryPtr->addChannel(1, 0.0, 0, 16, 13);
    slepEntryPtr->addChannel(1, 0.0, 0, -12, 11);
    slepEntryPtr->addChannel(1, 0.0, 0, -12, 13);
    slepEntryPtr->addChannel(1, 0.0, 0, -12, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, -14, 11);
    slepEntryPtr->addChannel(1, 0.0, 0, -14, 13);
    slepEntryPtr->addChannel(1, 0.0, 0, -14, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, -2, 1);
    slepEntryPtr->addChannel(1, 0.0, 0, -2, 3);
    slepEntryPtr->addChannel(1, 0.0, 0, -2, 5);
    slepEntryPtr->addChannel(1, 0.0, 0, -4, 1);
    slepEntryPtr->addChannel(1, 0.0, 0, -4, 3);
    slepEntryPtr->addChannel(1, 0.0, 0, -4, 5);
    slepEntryPtr->addChannel(1, 0.0, 0, -6, 1);
    slepEntryPtr->addChannel(1, 0.0, 0, -6, 3);
    slepEntryPtr->addChannel(1, 0.0, 0, -6, 5);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000022, 111, 16);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000022, 113, 16);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000022, 900111, 16);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000022, 16, 12, 11);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000022, 16, 14, 13);

   } else {
    slepEntryPtr->addChannel(1, 0.0, 0, 1000024, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000037, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000022, 16);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000023, 16);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000025, 16);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000035, 16);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000015, 24);
    slepEntryPtr->addChannel(1, 0.0, 0, 2000015, 24);
    slepEntryPtr->addChannel(1, 0.0, 0, 1000015, 37);
    slepEntryPtr->addChannel(1, 0.0, 0, 2000015, 37);
    slepEntryPtr->addChannel(1, 0.0, 0, -11, 11);
    slepEntryPtr->addChannel(1, 0.0, 0, -11, 13);
    slepEntryPtr->addChannel(1, 0.0, 0, -11, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, -13, 11);
    slepEntryPtr->addChannel(1, 0.0, 0, -13, 13);
    slepEntryPtr->addChannel(1, 0.0, 0, -13, 15);
    slepEntryPtr->addChannel(1, 0.0, 0, -1, 1);
    slepEntryPtr->addChannel(1, 0.0, 0, -1, 3);
    slepEntryPtr->addChannel(1, 0.0, 0, -1, 5);
    slepEntryPtr->addChannel(1, 0.0, 0, -3, 1);
    slepEntryPtr->addChannel(1, 0.0, 0, -3, 3);
    slepEntryPtr->addChannel(1, 0.0, 0, -3, 5);
    slepEntryPtr->addChannel(1, 0.0, 0, -5, 1);
    slepEntryPtr->addChannel(1, 0.0, 0, -5, 3);
    slepEntryPtr->addChannel(1, 0.0, 0, -5, 5);
  }

  return true;
}

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceSlepton::initConstants() {

  // Locally stored properties and couplings.
  s2W = coupSUSYPtr->sin2W;

  // Initialize functions for calculating 3-body widths
  stauWidths.setPointers(particleDataPtr,coupSUSYPtr,infoPtr);

}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceSlepton::calcPreFac(bool) {

  // Common coupling factors.
  alpEM  = coupSUSYPtr->alphaEM(mHat * mHat);
  preFac = 1.0 / (s2W * pow(mHat,3));

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceSlepton::calcWidth(bool) {

  // Slepton type -- in u_i/d_i and generation
  int ksusy = 1000000;
  int isl = (abs(idRes)/ksusy == 2) ? (abs(idRes)%10+1)/2 + 3
                                    : (abs(idRes)%10+1)/2;
  int il = (id2Abs-9)/2;
  bool islep = abs(idRes)%2;

  // Check that mass is above threshold.
  if (ps == 0.) return;
  widNow = 0.0;

  double fac = 0.0 , wid = 0.0;

  if (mult == 2) {
    // Two-body decays

    kinFac = (mHat * mHat - mf1 * mf1 - mf2 * mf2);
    fac = kinFac / (16.0 * M_PI * pow(mHat,3));

    // Case 1: RPV decays
    if (id1Abs < 17 && id2Abs < 17)  {

      wid = 0;
      int il2 = (id1Abs - 9)/2;

      //Case 1a: RPV LLE
      if(id1Abs > 10 && id2Abs > 10) {
        if (!coupSUSYPtr->isLLE) { widNow = 0.0; return;}

        if (!islep){ // sneutrino
          for (int isl2=1; isl2<3; isl2++)
            wid += norm(coupSUSYPtr->Rsv[isl][isl2]
                 * coupSUSYPtr->rvLLE[il][isl2][il2]);
        } else {
          for (int isl2=1; isl2<3; isl2++)
            wid += norm(coupSUSYPtr->Rsl[isl][isl2+3]
                 * coupSUSYPtr->rvLLE[isl2][il][il2]);
        }
      }
      //Case 1b: RPV LQD
      if(id1Abs < 10 && id2Abs < 10) {
        if (!coupSUSYPtr->isLQD) { widNow = 0.0; return;}
        if (!islep){ // sneutrino
          for (int isl2=1; isl2<3; isl2++)
            wid += norm(coupSUSYPtr->Rsv[isl][isl2]
                 * coupSUSYPtr->rvLQD[isl2][id1Abs][id2Abs]);
          wid *= 3.0; // colour factor

        } else {
          for (int isl2=1; isl2<3; isl2++)
            wid += norm(coupSUSYPtr->Rsl[isl][isl2+3]
                 * coupSUSYPtr->rvLLE[isl2][id1Abs][id2Abs]);
          wid *= 3.0; // colour factor
        }
      }
    }

    //Case 2: lepton + gaugino

    if (id1Abs > ksusy && id2Abs > 10 && id2Abs < 17) {
      for (int i=1; i<6 ; i++){
        // ~ell/~nu -> ~chi0 + ell/nu
        if (coupSUSYPtr->idNeut(i)==id1Abs && idRes%2 == id2Abs%2){
          fac = alpEM *  preFac / (2.0 * (1 - s2W));
          if (islep)
            wid = kinFac * (norm(coupSUSYPtr->LsllX[isl][il][i])
                + norm(coupSUSYPtr->RsllX[isl][il][i]))
                - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsllX[isl][il][i]
                * conj(coupSUSYPtr->RsllX[isl][il][i]));
          else
            wid = kinFac * (norm(coupSUSYPtr->LsvvX[isl][il][i])
                + norm(coupSUSYPtr->RsvvX[isl][il][i]))
                - 4.0 * mHat * mf2 * real(coupSUSYPtr->LsvvX[isl][il][i]
                * conj(coupSUSYPtr->RsvvX[isl][il][i]));
        }

        // ~ell/~nu -> ~chi- + nu/ell
        else if (i < 3 && coupSUSYPtr->idChar(i)==id1Abs
          && idRes%2 != id2Abs%2){

          fac = alpEM *  preFac / (4.0 * (1 - s2W));
          if (islep)
            wid = kinFac * (norm(coupSUSYPtr->LslvX[isl][il][i])
                + norm(coupSUSYPtr->RslvX[isl][il][i]))
                - 4.0 * mHat * mf2 * real(coupSUSYPtr->LslvX[isl][il][i]
                * conj(coupSUSYPtr->RslvX[isl][il][i]));
          else
            wid = kinFac * (norm(coupSUSYPtr->LslvX[isl][il][i])
                + norm(coupSUSYPtr->RslvX[isl][il][i]))
                - 4.0 * mHat * mf2 * real(coupSUSYPtr->LslvX[isl][il][i]
                * conj(coupSUSYPtr->RslvX[isl][il][i]));
        }
      }
    }

    //Case 3: ~l_i -> ~l_j + Z/W
    else if (id1Abs > ksusy+10 && id1Abs%100 < 17
      && (id2Abs == 23 || id2Abs == 24)){

      // factor of lambda^(3/2) = ps^3;
      fac = alpEM * preFac/(16.0 * pow2(mf2) * (1.0 - s2W)) * pow2(ps) ;

      int isl2 = (id1Abs/ksusy == 2) ? (id1Abs%10+1)/2 + 3: (id1Abs%10+1)/2;

      if (id2Abs == 23 && id1Abs%2 == idRes%2){
        if (islep)
          wid = norm(coupSUSYPtr->LslslZ[isl][isl2]
              + coupSUSYPtr->RslslZ[isl][isl2]);
        else
          wid = norm(coupSUSYPtr->LsvsvZ[isl][isl2]
              + coupSUSYPtr->RsvsvZ[isl][isl2]);
      }
      else if (id2Abs == 24 && id1Abs%2 != idRes%2){
        if (islep)
          wid = norm(coupSUSYPtr->LslsvW[isl2][isl]);
        else
          wid = norm(coupSUSYPtr->LslsvW[isl][isl2]);
      }
    }

    // TODO: Case ~l_i -> ~l_j + h/H

    widNow = fac * wid * ps * pow2(mHat);

    if (DBSUSY) cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs
                    <<" Width: "<<widNow<<endl;
  }
  else { // Case 4: special many-body stau decays
    // Case 4a: ~tau -> pi0 nu_tau ~chi0
    // Case 4b: ~tau -> rho/a1 nu_tau ~chi0
    // Case 4c: ~tau -> l nu_l nu_tau ~chi0

    // Check that there is a stau component
    double staufac = norm(coupSUSYPtr->Rsl[isl][3])
                   + norm(coupSUSYPtr->Rsl[isl][6]);
    if (staufac < 1.0e-6 || abs(mRes - particleDataPtr->m0(15)) > 0.0 ) return;
    if(id2Abs > 17)
      widNow = stauWidths.getWidth(idRes, id2Abs);
    else
      widNow = stauWidths.getWidth(idRes, id3Abs);

    widNow *= staufac;

  }

  return;

}

//==========================================================================

} //end namespace Pythia8
