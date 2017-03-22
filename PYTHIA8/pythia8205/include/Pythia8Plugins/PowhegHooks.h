// PowhegHooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Richard Corke, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Richard Corke.
// This class is used to perform a vetoed shower, where emissions
// already covered in a POWHEG NLO generator should be omitted.
// To first approximation the handover should happen at the SCALE
// of the LHA, but since the POWHEG-BOX uses a different pT definition
// than PYTHIA, both for ISR and FSR, a more sophisticated treatment
// is needed. See the online manual on POWHEG merging for details.

#ifndef Pythia8_PowhegHooks_H
#define Pythia8_PowhegHooks_H

// Includes
#include "Pythia8/Pythia.h"

namespace Pythia8 {

//==========================================================================

// Use userhooks to veto PYTHIA emissions above the POWHEG scale.

class PowhegHooks : public UserHooks {

public:

  // Constructor and destructor.
   PowhegHooks() {}
  ~PowhegHooks() {}

//--------------------------------------------------------------------------

  // Initialize settings, detailing merging strategy to use.
  bool initAfterBeams() {
    nFinal      = settingsPtr->mode("POWHEG:nFinal");
    vetoMode    = settingsPtr->mode("POWHEG:veto");
    vetoCount   = settingsPtr->mode("POWHEG:vetoCount");
    pThardMode  = settingsPtr->mode("POWHEG:pThard");
    pTemtMode   = settingsPtr->mode("POWHEG:pTemt");
    emittedMode = settingsPtr->mode("POWHEG:emitted");
    pTdefMode   = settingsPtr->mode("POWHEG:pTdef");
    MPIvetoMode = settingsPtr->mode("POWHEG:MPIveto");
    return true;
  }

//--------------------------------------------------------------------------

  // Routines to calculate the pT (according to pTdefMode) in a splitting:
  //   ISR: i (radiator after)  -> j (emitted after) k (radiator before)
  //   FSR: i (radiator before) -> j (emitted after) k (radiator after)
  // For the Pythia pT definition, a recoiler (after) must be specified.

  // Compute the Pythia pT separation. Based on pTLund function in History.cc
  inline double pTpythia(const Event &e, int RadAfterBranch, 
    int EmtAfterBranch, int RecAfterBranch, bool FSR) {

    // Convenient shorthands for later
    Vec4 radVec = e[RadAfterBranch].p();
    Vec4 emtVec = e[EmtAfterBranch].p();
    Vec4 recVec = e[RecAfterBranch].p();
    int  radID  = e[RadAfterBranch].id();

    // Calculate virtuality of splitting
    double sign = (FSR) ? 1. : -1.;
    Vec4 Q(radVec + sign * emtVec);
    double Qsq = sign * Q.m2Calc();

    // Mass term of radiator
    double m2Rad = (abs(radID) >= 4 && abs(radID) < 7) ?
                   pow2(particleDataPtr->m0(radID)) : 0.;

    // z values for FSR and ISR
    double z, pTnow;
    if (FSR) {
      // Construct 2 -> 3 variables
      Vec4 sum = radVec + recVec + emtVec;
      double m2Dip = sum.m2Calc();
      double x1 = 2. * (sum * radVec) / m2Dip;
      double x3 = 2. * (sum * emtVec) / m2Dip;
      z     = x1 / (x1 + x3);
      pTnow = z * (1. - z);

    } else {
      // Construct dipoles before/after splitting
      Vec4 qBR(radVec - emtVec + recVec);
      Vec4 qAR(radVec + recVec);
      z     = qBR.m2Calc() / qAR.m2Calc();
      pTnow = (1. - z);
    }

    // Virtuality with correct sign
    pTnow *= (Qsq - sign * m2Rad);

    // Can get negative pT for massive splittings
    if (pTnow < 0.) {
      cout << "Warning: pTpythia was negative" << endl;
      return -1.;
    }

    // Return pT
    return sqrt(pTnow);
  }

  // Compute the POWHEG pT separation between i and j
  inline double pTpowheg(const Event &e, int i, int j, bool FSR) {

    // pT value for FSR and ISR
    double pTnow = 0.;
    if (FSR) {
      // POWHEG d_ij (in CM frame). Note that the incoming beams have not
      // been updated in the parton systems pointer yet (i.e. prior to any
      // potential recoil).
      int iInA = partonSystemsPtr->getInA(0);
      int iInB = partonSystemsPtr->getInB(0);
      double betaZ = - ( e[iInA].pz() + e[iInB].pz() ) /
                       ( e[iInA].e()  + e[iInB].e()  );
      Vec4 iVecBst(e[i].p()), jVecBst(e[j].p());
      iVecBst.bst(0., 0., betaZ);
      jVecBst.bst(0., 0., betaZ);
      pTnow = sqrt( (iVecBst + jVecBst).m2Calc() *
                    iVecBst.e() * jVecBst.e() /
                    pow2(iVecBst.e() + jVecBst.e()) );

    } else {
      // POWHEG pT_ISR is just kinematic pT
      pTnow = e[j].pT();
    }

    // Check result
    if (pTnow < 0.) {
      cout << "Warning: pTpowheg was negative" << endl;
      return -1.;
    }

    return pTnow;
  }

  // Calculate pT for a splitting based on pTdefMode.
  // If j is -1, all final-state partons are tried.
  // If i, k, r and xSR are -1, then all incoming and outgoing
  // partons are tried.
  // xSR set to 0 means ISR, while xSR set to 1 means FSR
  inline double pTcalc(const Event &e, int i, int j, int k, int r, int xSRin) {

    // Loop over ISR and FSR if necessary
    double pTemt = -1., pTnow;
    int xSR1 = (xSRin == -1) ? 0 : xSRin;
    int xSR2 = (xSRin == -1) ? 2 : xSRin + 1;
    for (int xSR = xSR1; xSR < xSR2; xSR++) {
      // FSR flag
      bool FSR = (xSR == 0) ? false : true;

      // If all necessary arguments have been given, then directly calculate.
      // POWHEG ISR and FSR, need i and j.
      if ((pTdefMode == 0 || pTdefMode == 1) && i > 0 && j > 0) {
        pTemt = pTpowheg(e, i, j, (pTdefMode == 0) ? false : FSR);

      // Pythia ISR, need i, j and r.
      } else if (!FSR && pTdefMode == 2 && i > 0 && j > 0 && r > 0) {
        pTemt = pTpythia(e, i, j, r, FSR);

      // Pythia FSR, need k, j and r.
      } else if (FSR && pTdefMode == 2 && j > 0 && k > 0 && r > 0) {
        pTemt = pTpythia(e, k, j, r, FSR);

      // Otherwise need to try all possible combinations.
      } else {
        // Start by finding incoming legs to the hard system after
        // branching (radiator after branching, i for ISR).
        // Use partonSystemsPtr to find incoming just prior to the
        // branching and track mothers.
        int iInA = partonSystemsPtr->getInA(0);
        int iInB = partonSystemsPtr->getInB(0);
        while (e[iInA].mother1() != 1) { iInA = e[iInA].mother1(); }
        while (e[iInB].mother1() != 2) { iInB = e[iInB].mother1(); }

        // If we do not have j, then try all final-state partons
        int jNow = (j > 0) ? j : 0;
        int jMax = (j > 0) ? j + 1 : e.size();
        for (; jNow < jMax; jNow++) {

          // Final-state and coloured jNow only
          if (!e[jNow].isFinal() || e[jNow].colType() == 0) continue;

          // POWHEG
          if (pTdefMode == 0 || pTdefMode == 1) {

            // ISR - only done once as just kinematical pT
            if (!FSR) {
              pTnow = pTpowheg(e, iInA, jNow, (pTdefMode == 0) ? false : FSR);
              if (pTnow > 0.) pTemt = (pTemt < 0) ? pTnow : min(pTemt, pTnow);

              // FSR - try all outgoing partons from system before branching
              // as i. Note that for the hard system, there is no
              // "before branching" information.
              } else {

                int outSize = partonSystemsPtr->sizeOut(0);
                for (int iMem = 0; iMem < outSize; iMem++) {
                  int iNow = partonSystemsPtr->getOut(0, iMem);

                  // Coloured only, i != jNow and no carbon copies
                  if (iNow == jNow || e[iNow].colType() == 0) continue;
                  if (jNow == e[iNow].daughter1()
                    && jNow == e[iNow].daughter2()) continue;

                  pTnow = pTpowheg(e, iNow, jNow, (pTdefMode == 0)
                    ? false : FSR);
                  if (pTnow > 0.) pTemt = (pTemt < 0)
                    ? pTnow : min(pTemt, pTnow);
                } // for (iMem)

              } // if (!FSR)

          // Pythia
          } else if (pTdefMode == 2) {

            // ISR - other incoming as recoiler
            if (!FSR) {
              pTnow = pTpythia(e, iInA, jNow, iInB, FSR);
              if (pTnow > 0.) pTemt = (pTemt < 0) ? pTnow : min(pTemt, pTnow);
              pTnow = pTpythia(e, iInB, jNow, iInA, FSR);
              if (pTnow > 0.) pTemt = (pTemt < 0) ? pTnow : min(pTemt, pTnow);

            // FSR - try all final-state coloured partons as radiator
            //       after emission (k).
            } else {
              for (int kNow = 0; kNow < e.size(); kNow++) {
                if (kNow == jNow || !e[kNow].isFinal() ||
                    e[kNow].colType() == 0) continue;

                // For this kNow, need to have a recoiler.
                // Try two incoming.
                pTnow = pTpythia(e, kNow, jNow, iInA, FSR);
                if (pTnow > 0.) pTemt = (pTemt < 0)
                  ? pTnow : min(pTemt, pTnow);
                pTnow = pTpythia(e, kNow, jNow, iInB, FSR);
                if (pTnow > 0.) pTemt = (pTemt < 0)
                  ? pTnow : min(pTemt, pTnow);

                // Try all other outgoing.
                for (int rNow = 0; rNow < e.size(); rNow++) {
                  if (rNow == kNow || rNow == jNow ||
                      !e[rNow].isFinal() || e[rNow].colType() == 0) continue;
                  pTnow = pTpythia(e, kNow, jNow, rNow, FSR);
                  if (pTnow > 0.) pTemt = (pTemt < 0)
                    ? pTnow : min(pTemt, pTnow);
                } // for (rNow)

              } // for (kNow)
            } // if (!FSR)
          } // if (pTdefMode)
        } // for (j)
      }
    } // for (xSR)

    return pTemt;
  }

//--------------------------------------------------------------------------

  // Extraction of pThard based on the incoming event.
  // Assume that all the final-state particles are in a continuous block
  // at the end of the event and the final entry is the POWHEG emission.
  // If there is no POWHEG emission, then pThard is set to SCALUP.

  inline bool canVetoMPIStep()    { return true; }
  inline int  numberVetoMPIStep() { return 1; }
  inline bool doVetoMPIStep(int nMPI, const Event &e) {
    // Extra check on nMPI
    if (nMPI > 1) return false;

    // Find if there is a POWHEG emission. Go backwards through the
    // event record until there is a non-final particle. Also sum pT and
    // find pT_1 for possible MPI vetoing
    int    count = 0;
    double pT1 = 0., pTsum = 0.;
    for (int i = e.size() - 1; i > 0; i--) {
      if (e[i].isFinal()) {
        count++;
        pT1    = e[i].pT();
        pTsum += e[i].pT();
      } else break;
    }
    // Extra check that we have the correct final state
    if (count != nFinal && count != nFinal + 1) {
      cout << "Error: wrong number of final state particles in event" << endl;
      exit(1);
    }
    // Flag if POWHEG radiation present and index
    bool isEmt = (count == nFinal) ? false : true;
    int  iEmt  = (isEmt) ? e.size() - 1 : -1;

    // If there is no radiation or if pThardMode is 0 then set pThard = SCALUP.
    if (!isEmt || pThardMode == 0) {
      pThard = infoPtr->scalup();

    // If pThardMode is 1 then the pT of the POWHEG emission is checked against
    // all other incoming and outgoing partons, with the minimal value taken
    } else if (pThardMode == 1) {
      pThard = pTcalc(e, -1, iEmt, -1, -1, -1);

    // If pThardMode is 2, then the pT of all final-state partons is checked
    // against all other incoming and outgoing partons, with the minimal value
    // taken
    } else if (pThardMode == 2) {
      pThard = pTcalc(e, -1, -1, -1, -1, -1);
    }

    // Find MPI veto pT if necessary
    if (MPIvetoMode == 1) {
      pTMPI = (isEmt) ? pTsum / 2. : pT1;
    }

    // Initialise other variables
    accepted   = false;
    nAcceptSeq = nISRveto = nFSRveto = 0;

    // Do not veto the event
    return false;
  }

//--------------------------------------------------------------------------

  // ISR veto

  inline bool canVetoISREmission() { return (vetoMode == 0) ? false : true; }
  inline bool doVetoISREmission(int, const Event &e, int iSys) {
    // Must be radiation from the hard system
    if (iSys != 0) return false;

    // If we already have accepted 'vetoCount' emissions in a row, do nothing
    if (vetoMode == 1 && nAcceptSeq >= vetoCount) return false;

    // Pythia radiator after, emitted and recoiler after.
    int iRadAft = -1, iEmt = -1, iRecAft = -1;
    for (int i = e.size() - 1; i > 0; i--) {
      if      (iRadAft == -1 && e[i].status() == -41) iRadAft = i;
      else if (iEmt    == -1 && e[i].status() ==  43) iEmt    = i;
      else if (iRecAft == -1 && e[i].status() == -42) iRecAft = i;
      if (iRadAft != -1 && iEmt != -1 && iRecAft != -1) break;
    }
    if (iRadAft == -1 || iEmt == -1 || iRecAft == -1) {
      e.list();
      cout << "Error: couldn't find Pythia ISR emission" << endl;
      exit(1);
    }

    // pTemtMode == 0: pT of emitted w.r.t. radiator
    // pTemtMode == 1: min(pT of emitted w.r.t. all incoming/outgoing)
    // pTemtMode == 2: min(pT of all outgoing w.r.t. all incoming/outgoing)
    int xSR      = (pTemtMode == 0) ? 0       : -1;
    int i        = (pTemtMode == 0) ? iRadAft : -1;
    int j        = (pTemtMode != 2) ? iEmt    : -1;
    int k        = -1;
    int r        = (pTemtMode == 0) ? iRecAft : -1;
    double pTemt = pTcalc(e, i, j, k, r, xSR);

    // Veto if pTemt > pThard
    if (pTemt > pThard) {
      nAcceptSeq = 0;
      nISRveto++;
      return true;
    }

    // Else mark that an emission has been accepted and continue
    nAcceptSeq++;
    accepted = true;
    return false;
  }

//--------------------------------------------------------------------------

  // FSR veto

  inline bool canVetoFSREmission() { return (vetoMode == 0) ? false : true; }
  inline bool doVetoFSREmission(int, const Event &e, int iSys, bool) {
    // Must be radiation from the hard system
    if (iSys != 0) return false;

    // If we already have accepted 'vetoCount' emissions in a row, do nothing
    if (vetoMode == 1 && nAcceptSeq >= vetoCount) return false;

    // Pythia radiator (before and after), emitted and recoiler (after)
    int iRecAft = e.size() - 1;
    int iEmt    = e.size() - 2;
    int iRadAft = e.size() - 3;
    int iRadBef = e[iEmt].mother1();
    if ( (e[iRecAft].status() != 52 && e[iRecAft].status() != -53) ||
         e[iEmt].status() != 51 || e[iRadAft].status() != 51) {
      e.list();
      cout << "Error: couldn't find Pythia FSR emission" << endl;
      exit(1);
    }

    // Behaviour based on pTemtMode:
    //  0 - pT of emitted w.r.t. radiator before
    //  1 - min(pT of emitted w.r.t. all incoming/outgoing)
    //  2 - min(pT of all outgoing w.r.t. all incoming/outgoing)
    int xSR = (pTemtMode == 0) ? 1       : -1;
    int i   = (pTemtMode == 0) ? iRadBef : -1;
    int k   = (pTemtMode == 0) ? iRadAft : -1;
    int r   = (pTemtMode == 0) ? iRecAft : -1;

    // When pTemtMode is 0 or 1, iEmt has been selected
    double pTemt = 0.;
    if (pTemtMode == 0 || pTemtMode == 1) {
      // Which parton is emitted, based on emittedMode:
      //  0 - Pythia definition of emitted
      //  1 - Pythia definition of radiated after emission
      //  2 - Random selection of emitted or radiated after emission
      //  3 - Try both emitted and radiated after emission
      int j = iRadAft;
      if (emittedMode == 0 || (emittedMode == 2 && rndmPtr->flat() < 0.5)) j++;

      for (int jLoop = 0; jLoop < 2; jLoop++) {
        if      (jLoop == 0) pTemt = pTcalc(e, i, j, k, r, xSR);
        else if (jLoop == 1) pTemt = min(pTemt, pTcalc(e, i, j, k, r, xSR));

        // For emittedMode == 3, have tried iRadAft, now try iEmt
        if (emittedMode != 3) break;
        if (k != -1) swap(j, k); else j = iEmt;
      }

    // If pTemtMode is 2, then try all final-state partons as emitted
    } else if (pTemtMode == 2) {
      pTemt = pTcalc(e, i, -1, k, r, xSR);

    }

    // Veto if pTemt > pThard
    if (pTemt > pThard) {
      nAcceptSeq = 0;
      nFSRveto++;
      return true;
    }

    // Else mark that an emission has been accepted and continue
    nAcceptSeq++;
    accepted = true;
    return false;
  }

//--------------------------------------------------------------------------

  // MPI veto

  inline bool canVetoMPIEmission() {return (MPIvetoMode == 0) ? false : true;}
  inline bool doVetoMPIEmission(int, const Event &e) {
    if (MPIvetoMode == 1) {
      if (e[e.size() - 1].pT() > pTMPI) return true;
    }
    return false;
  }

//--------------------------------------------------------------------------

  // Functions to return information

  inline int    getNISRveto() { return nISRveto; }
  inline int    getNFSRveto() { return nFSRveto; }

//--------------------------------------------------------------------------

  // Private members

private:
  int    nFinal, vetoMode, vetoCount, pThardMode, pTemtMode,
         emittedMode, pTdefMode, MPIvetoMode;
  double pThard, pTMPI;
  bool   accepted;
  // The number of accepted emissions (in a row)
  int nAcceptSeq;
  // Statistics on vetos
  unsigned long int nISRveto, nFSRveto;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_PowhegHooks_H

