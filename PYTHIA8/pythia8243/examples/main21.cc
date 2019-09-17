// main21.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how to feed in a single particle (including a resonance)
// or a toy parton-level configurations.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Single-particle gun. The particle must be a colour singlet.
// Input: flavour, energy, direction (theta, phi).
// If theta < 0 then random choice over solid angle.
// Optional final argument to put particle at rest => E = m.

void fillParticle(int id, double ee, double thetaIn, double phiIn,
  Event& event, ParticleData& pdt, Rndm& rndm, bool atRest = false,
  bool hasLifetime = false) {

  // Reset event record to allow for new event.
  event.reset();

  // Select particle mass; where relevant according to Breit-Wigner.
  double mm = pdt.mSel(id);
  double pp = sqrtpos(ee*ee - mm*mm);

  // Special case when particle is supposed to be at rest.
  if (atRest) {
    ee = mm;
    pp = 0.;
  }

  // Angles as input or uniform in solid angle.
  double cThe, sThe, phi;
  if (thetaIn >= 0.) {
    cThe = cos(thetaIn);
    sThe = sin(thetaIn);
    phi  = phiIn;
  } else {
    cThe = 2. * rndm.flat() - 1.;
    sThe = sqrtpos(1. - cThe * cThe);
    phi = 2. * M_PI * rndm.flat();
  }

  // Store the particle in the event record.
  int iNew = event.append( id, 1, 0, 0, pp * sThe * cos(phi),
    pp * sThe * sin(phi), pp * cThe, ee, mm);

  // Generate lifetime, to give decay away from primary vertex.
  if (hasLifetime) event[iNew].tau( event[iNew].tau0() * rndm.exp() );

}

//==========================================================================

// Simple method to do the filling of partons into the event record.

void fillPartons(int type, double ee, Event& event, ParticleData& pdt,
  Rndm& rndm) {

  // Reset event record to allow for new event.
  event.reset();

  // Information on a q qbar system, to be hadronized.
  if (type == 1 || type == 12) {
    int    id = 2;
    double mm = pdt.m0(id);
    double pp = sqrtpos(ee*ee - mm*mm);
    event.append(  id, 23, 101,   0, 0., 0.,  pp, ee, mm);
    event.append( -id, 23,   0, 101, 0., 0., -pp, ee, mm);

  // Information on a g g system, to be hadronized.
  } else if (type == 2 || type == 13) {
    event.append( 21, 23, 101, 102, 0., 0.,  ee, ee);
    event.append( 21, 23, 102, 101, 0., 0., -ee, ee);

  // Information on a g g g system, to be hadronized.
  } else if (type == 3) {
    event.append( 21, 23, 101, 102,        0., 0.,        ee, ee);
    event.append( 21, 23, 102, 103,  0.8 * ee, 0., -0.6 * ee, ee);
    event.append( 21, 23, 103, 101, -0.8 * ee, 0., -0.6 * ee, ee);

  // Information on a q q q junction system, to be hadronized.
  } else if (type == 4 || type == 5) {

    // Need a colour singlet mother parton to define junction origin.
    event.append( 1000022, -21, 0, 0, 2, 4, 0, 0,
                  0., 0., 1.01 * ee, 1.01 * ee);

    // The three endpoint q q q; the minimal system.
    double rt75 = sqrt(0.75);
    event.append( 2, 23, 1, 0, 0, 0, 101, 0,
                          0., 0., 1.01 * ee, 1.01 * ee);
    event.append( 2, 23, 1, 0, 0, 0, 102, 0,
                   rt75 * ee, 0., -0.5 * ee,        ee );
    event.append( 1, 23, 1, 0, 0, 0, 103, 0,
                  -rt75 * ee, 0., -0.5 * ee,        ee );

    // Define the qqq configuration as starting point for adding gluons.
    if (type == 5) {
      int colNow[4] = {0, 101, 102, 103};
      Vec4 pQ[4];
      pQ[1] = Vec4(0., 0., 1., 0.);
      pQ[2] = Vec4( rt75, 0., -0.5, 0.);
      pQ[3] = Vec4(-rt75, 0., -0.5, 0.);

      // Minimal cos(q-g opening angle), allows more or less nasty events.
      double cosThetaMin =0.;

      // Add a few gluons (almost) at random.
      for (int nglu = 0; nglu < 5; ++nglu) {
        int iq = 1 + int( 2.99999 * rndm.flat() );
        double px, py, pz, e, prod;
        do {
          e =  ee * rndm.flat();
          double cThe = 2. * rndm.flat() - 1.;
          double phi = 2. * M_PI * rndm.flat();
          px = e * sqrt(1. - cThe*cThe) * cos(phi);
          py = e * sqrt(1. - cThe*cThe) * sin(phi);
          pz = e * cThe;
          prod = ( px * pQ[iq].px() + py * pQ[iq].py() + pz * pQ[iq].pz() )
            / e;
        } while (prod < cosThetaMin);
        int colNew = 104 + nglu;
        event.append( 21, 23, 1, 0, 0, 0, colNew, colNow[iq],
          px, py, pz, e, 0.);
        colNow[iq] = colNew;
      }
      // Update daughter range of mother.
      event[1].daughters(2, event.size() - 1);

    }

  // Information on a q q qbar qbar dijunction system, to be hadronized.
  } else if (type >= 6) {

    // The two fictitious beam remnant particles; needed for junctions.
    event.append( 2212, -12, 0, 0, 3, 5, 0, 0, 0., 0., ee, ee, 0.);
    event.append(-2212, -12, 0, 0, 6, 8, 0, 0, 0., 0., ee, ee, 0.);

    // Opening angle between "diquark" legs.
    double theta = 0.2;
    double cThe = cos(theta);
    double sThe = sin(theta);

    // Set one colour depending on whether more gluons or not.
    int acol = (type == 6) ? 103 : 106;

    // The four endpoint q q qbar qbar; the minimal system.
    // Two additional fictitious partons to make up original beams.
    event.append(  2,   23, 1, 0, 0, 0, 101, 0,
                  ee * sThe, 0.,  ee * cThe, ee, 0.);
    event.append(  1,   23, 1, 0, 0, 0, 102, 0,
                 -ee * sThe, 0.,  ee * cThe, ee, 0.);
    event.append(  2, -21, 1, 0, 0, 0, 103, 0,
                         0., 0.,  ee       , ee, 0.);
    event.append( -2,   23, 2, 0, 0, 0, 0, 104,
                  ee * sThe, 0., -ee * cThe, ee, 0.);
    event.append( -1,   23, 2, 0, 0, 0, 0, 105,
                 -ee * sThe, 0., -ee * cThe, ee, 0.);
    event.append( -2, -21, 2, 0, 0, 0, 0, acol,
                         0., 0., -ee       , ee, 0.);

    // Add extra gluons on string between junctions.
    if (type == 7) {
      event.append( 21, 23, 5, 8, 0, 0, 103, 106, 0., ee, 0., ee, 0.);
    } else if (type == 8) {
      event.append( 21, 23, 5, 8, 0, 0, 103, 108, 0., ee, 0., ee, 0.);
      event.append( 21, 23, 5, 8, 0, 0, 108, 106, 0.,-ee, 0., ee, 0.);
    } else if (type == 9) {
      event.append( 21, 23, 5, 8, 0, 0, 103, 107, 0., ee, 0., ee, 0.);
      event.append( 21, 23, 5, 8, 0, 0, 107, 108, ee, 0., 0., ee, 0.);
      event.append( 21, 23, 5, 8, 0, 0, 108, 106, 0.,-ee, 0., ee, 0.);
    } else if (type == 10) {
      event.append( 21, 23, 5, 8, 0, 0, 103, 107, 0., ee, 0., ee, 0.);
      event.append( 21, 23, 5, 8, 0, 0, 107, 108, ee, 0., 0., ee, 0.);
      event.append( 21, 23, 5, 8, 0, 0, 108, 109, 0.,-ee, 0., ee, 0.);
      event.append( 21, 23, 5, 8, 0, 0, 109, 106,-ee, 0., 0., ee, 0.);
    }

  // No more cases: done.
  }
}

//==========================================================================

int main() {

  // Pick kind of events to generate:
  // 0 = single-particle gun.
  // 1 = q qbar.
  // 2 = g g.
  // 3 = g g g.
  // 4 = minimal q q q junction topology.
  // 5 = q q q junction topology with gluons on the strings.
  // 6 = q q qbar qbar dijunction topology, no gluons.
  // 7 - 10 = ditto, but with 1 - 4 gluons on string between junctions.
  // 11 = single-resonance gun.
  // 12 = q qbar plus parton shower.
  // 13 = g g plus parton shower.
  int type = 11;

  // Set particle species and energy for single-particle gun.
  int    idGun  = (type == 0) ? 15 : 25;
  double eeGun  = (type == 0) ? 20. : 125.;
  bool   atRest = (type == 0) ? false : true;

  // The single-particle gun produces a particle at the origin, and by default
  // decays it there. When hasLifetime = true instead a finite lifetime is
  // selected and used to generate a displaced  decay vertex.
  bool   hasLifetime = (type == 0) ? true : false;

  // Set typical energy per parton.
  double ee = 20.0;

  // Set number of events to generate and to list.
  int nEvent = 10000;
  int nList = 3;

  // Generator; shorthand for event and particleData.
  Pythia pythia;
  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;

  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  pythia.readString("ProcessLevel:all = off");

  // Optionally switch off resonance decays, or only showers in them.
  //pythia.readString("ProcessLevel:resonanceDecays = off");
  //pythia.readString("PartonLevel:FSRinResonances = off");

  // Optionally switch off ordinary decays.
  //pythia.readString("HadronLevel:Decay = off");

  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Set true to also see space-time information in event listings.
  bool showScaleAndVertex = (type == 0) ? true : false;

  // Initialize.
  pythia.init();

  // Book histograms.
  Hist epCons("deviation from energy-momentum conservation",100,0.,1e-4);
  Hist chgCons("deviation from charge conservation",57,-9.5,9.5);
  Hist nFinal("final particle multiplicity",100,-0.5,99.5);
  Hist dnparticledp("dn/dp for particles",100,0.,ee);
  Hist status85("multiplicity status code 85",50,-0.5,49.5);
  Hist status86("multiplicity status code 86",50,-0.5,49.5);
  Hist status83("multiplicity status code 83",50,-0.5,49.5);
  Hist status84("multiplicity status code 84",50,-0.5,49.5);
  Hist dndtheta("particle flow in event plane",100,-M_PI,M_PI);
  Hist dedtheta("energy flow in event plane",100,-M_PI,M_PI);
  Hist dpartondtheta("parton flow in event plane",100,-M_PI,M_PI);
  Hist dndyAnti("dn/dy primaries antijunction",100, -10., 10.);
  Hist dndyJun("dn/dy primaries junction",100, -10., 10.);
  Hist dndySum("dn/dy all primaries",100, -10., 10.);

  // Begin of event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Set up single particle, with random direction in solid angle.
    if (type == 0 || type == 11) fillParticle( idGun, eeGun, -1., 0.,
      event, pdt, pythia.rndm, atRest, hasLifetime);

    // Set up parton-level configuration.
    else fillPartons( type, ee, event, pdt, pythia.rndm);

    // To have partons shower they must be set maximum allowed scale.
    // (Can be set individually to restrict radiation differently.)
    if (type == 12 || type == 13) {
      double scale = ee;
      event[1].scale( scale);
      event[2].scale( scale);

      // Now actually do the shower, for range of partons, and max scale.
      // (Most restrictive of global and individual applied to each parton.)
      pythia.forceTimeShower( 1, 2, scale);
    }

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // List first few events.
    if (iEvent < nList) {
      event.list(showScaleAndVertex);
      // Also list junctions.
      event.listJunctions();
    }

    // Initialize statistics.
    Vec4 pSum = - event[0].p();
    double chargeSum = 0.;
    if (type == 0) chargeSum = -event[1].charge();
    if (type == 4 || type == 5) chargeSum = -1;
    int nFin = 0;
    int n85 = 0;
    int n86 = 0;
    int n83 = 0;
    int n84 = 0;

    // Loop over all particles.
    for (int i = 0; i < event.size(); ++i) {
      int status = event[i].statusAbs();

      // Find any unrecognized particle codes.
      int id = event[i].id();
      if (id == 0 || !pdt.isParticle(id))
        cout << " Error! Unknown code id = " << id << "\n";

      // Find particles with E-p mismatch.
      double eCalc = event[i].eCalc();
      if (abs(eCalc/event[i].e() - 1.) > 1e-6) cout << " e mismatch, i = "
        << i << " e_nominal = " << event[i].e() << " e-from-p = "
        << eCalc << " m-from-e " << event[i].mCalc() << "\n";

      // Parton flow in event plane.
      if (status == 71 || status == 72) {
        double thetaXZ = event[i].thetaXZ();
        dpartondtheta.fill(thetaXZ);
      }

      // Origin of primary hadrons.
      if (status == 85) ++n85;
      if (status == 86) ++n86;
      if (status == 83) ++n83;
      if (status == 84) ++n84;

      // Flow of primary hadrons in the event plane.
      if (status > 80 && status < 90) {
        double eAbs = event[i].e();
        if (eAbs < 0.) {cout << " e < 0 line " << i; event.list();}
        double thetaXZ = event[i].thetaXZ();
        dndtheta.fill(thetaXZ);
        dedtheta.fill(thetaXZ, eAbs);

        // Rapidity distribution of primary hadrons.
        double y = event[i].y();
        dndySum.fill(y);
        if (type >= 6) {
          int motherId = event[event[i].mother1()].id();
          if (motherId > 0 ) dndyJun.fill(event[i].y());
          else dndyAnti.fill(event[i].y());
        }
      }

      // Study final-state particles.
      if (event[i].isFinal()) {
        pSum += event[i].p();
        chargeSum += event[i].charge();
        nFin++;
        double pAbs = event[i].pAbs();
        dnparticledp.fill(pAbs);
      }
    }

    // Fill histograms once for each event.
    double epDev = abs(pSum.e()) + abs(pSum.px()) + abs(pSum.py())
      + abs(pSum.pz());
    epCons.fill(epDev);
    chgCons.fill(chargeSum);
    nFinal.fill(nFin);
    status85.fill(n85);
    status86.fill(n86);
    status83.fill(n83);
    status84.fill(n84);
    if (epDev > 1e-3  || abs(chargeSum) > 0.1) event.list();

  // End of event loop.
  }

  // Print statistics, histograms and done.
  pythia.stat();
  cout << epCons << chgCons << nFinal << dnparticledp
       << dndtheta << dedtheta << dpartondtheta << dndySum;
  if (type >= 4) cout << status85 << status86 << status83
       << status84;
  if (type >= 6) cout << dndyJun << dndyAnti;

  // Done.
  return 0;
}
