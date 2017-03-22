// main18.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how to write an event filter.
// No new functionality is involved - all could be done in the main program
// - but the division of tasks may be more convenient for recurrent cuts.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// The EventFilter class.

// The constructor takes the following arguments
// select = 1 : keep only final particles.
//        = 2 : keep only final visible particles (i.e. not neutrinos).
//        = 3 : keep only final charged particles.
// etaMax (default = 50) : keep only particles with pseudorapidity
//        |eta| < etaMax.
// pTminCharged (default = 0) : keep a charged particle only if
//        its transverse momentum pT < pTminCharged.
// pTminNeutral (default = 0) : keep a neutral particle only if
//        its transverse momentum pT < pTminNeutral.

// Main methods:
// filter( event) takes an event record as input and analyzes it.
// size() returns the number of particles kept.
// index(i) returns the index in the full event of the i'th kept particle.
// particlePtr(i) returns a pointer to the i'th kept particle.
// particleRef(i) returns a reference to the i'th kept particle.
// list() gives a listing of the kept particles only.

class EventFilter {

public:

  // Constructor sets properties of filter.
  EventFilter( int selectIn, double etaMaxIn = 50.,
    double pTminChargedIn = 0., double pTminNeutralIn = 0.)
    : select(selectIn), etaMax(etaMaxIn), pTminCharged(pTminChargedIn),
    pTminNeutral(pTminNeutralIn) {}

  // Analysis of each new event to find acceptable particles.
  void filter(Event& event);

  // Return size of array, and index of a particle.
  int size()       const {return keptPtrs.size();}
  int index(int i) const {return keptIndx[i];}

  // Return pointer or reference to a particle.
  Particle* particlePtr(int i) {return  keptPtrs[i];}
  Particle& particleRef(int i) {return *keptPtrs[i];}

  // List kept particles only.
  void list(ostream& os = cout);

private:

  // Filter properties, set by constructor.
  int    select;
  double etaMax, pTminCharged, pTminNeutral;

  // Kept particle indices and pointers, referring to original event.
  vector<int>       keptIndx;
  vector<Particle*> keptPtrs;

};

//--------------------------------------------------------------------------

// The filter method.

void EventFilter::filter(Event& event) {

  // Reset arrays in preparation for new event.
  keptIndx.resize(0);
  keptPtrs.resize(0);

  // Loop over all particles in the event record.
  for (int i = 0; i < event.size(); ++i) {

    // Skip if particle kind selection criteria not fulfilled.
    if (!event[i].isFinal()) continue;
    if (select == 2 && !event[i].isVisible()) continue;
    bool isCharged = event[i].isCharged();
    if (select == 3 && !isCharged) continue;

    // Skip if too large pseudorapidity.
    if (abs(event[i].eta()) > etaMax) continue;

    // Skip if too small pT.
    if       (isCharged && event[i].pT() < pTminCharged) continue;
    else if (!isCharged && event[i].pT() < pTminNeutral) continue;

    // Add particle to vectors of indices and pointers.
    keptIndx.push_back( i );
    keptPtrs.push_back( &event[i] );

  // End of particle loop. Done.
  }

}

//--------------------------------------------------------------------------

// The list method: downscaled version of Event::list.

void EventFilter::list(ostream& os) {

  // Header.
  os << "\n --------  PYTHIA Event Listing  (filtered)  ------------------"
     << "-----------------------------------------------------------------"
     << "----\n \n    no        id   name            status     mothers  "
     << " daughters     colours      p_x        p_y        p_z         e  "
     << "        m \n";

  // At high energy switch to scientific format for momenta.
  double eSum = 0.;
  for (int iKept = 0; iKept < size(); ++iKept) eSum += keptPtrs[iKept]->e();
  bool useFixed = (eSum < 1e5);

  // Listing of kept particles in event.
  for (int iKept = 0; iKept < size(); ++iKept) {
    int i = keptIndx[iKept];
    Particle& pt = *keptPtrs[iKept];

    // Basic line for a particle, always printed.
    os << setw(6) << i << setw(10) << pt.id() << "   " << left
       << setw(18) << pt.nameWithStatus(18) << right << setw(4)
       << pt.status() << setw(6) << pt.mother1() << setw(6)
       << pt.mother2() << setw(6) << pt.daughter1() << setw(6)
       << pt.daughter2() << setw(6) << pt.col() << setw(6) << pt.acol()
       << ( (useFixed) ? fixed : scientific ) << setprecision(3)
       << setw(11) << pt.px() << setw(11) << pt.py() << setw(11)
       << pt.pz() << setw(11) << pt.e() << setw(11) << pt.m() << "\n";
  }

  // Listing finished.
  os << "\n --------  End PYTHIA Event Listing  ----------------------------"
     << "-------------------------------------------------------------------"
     << endl;
}


//==========================================================================

// Use the EventFilter method to plot some event properties.

int main() {

  // Number of events to generate, to list, to allow aborts.
  int    nEvent   = 100;
  int    nList    = 1;
  int    nAbort   = 3;

  // Declare generator.
  Pythia pythia;

  // Hard QCD events with pThat > 100.
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 100.");

  // No automatic event listings - do it manually below.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Initialization for LHC.
  pythia.init();

  // Values for filter.
  int    select   = 3;
  double etaMax   = 3.;
  double pTminChg = 1.;

  // Declare Event Filter according to specification.
  EventFilter filter( select, etaMax, pTminChg);

  // Histograms.
  Hist nCharged(   "selected charged multiplicity",     100, -0.5, 199.5);
  Hist etaCharged( "selected charged eta distribution", 100, -5.0, 5.0);
  Hist pTCharged(  "selected charged pT distribution",  100,  0.0, 50.0);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Find final charged particles with |eta| < 3 and pT > 1 GeV.
    filter.filter( pythia.event);

    // List first few events, both complete and after filtering.
    if (iEvent < nList) {
      pythia.info.list();
      pythia.process.list();
      pythia.event.list();
      filter.list();
    }

    // Analyze selected particle sample.
    nCharged.fill( filter.size() );
    for (int i = 0; i < filter.size(); ++i) {
      // Use both reference and pointer notation to illustrate freedom.
      etaCharged.fill( filter.particleRef(i).eta() );
      pTCharged.fill(  filter.particlePtr(i)->pT() );
    }

  // End of event loop.
  }

  // Final statistics.
  pythia.stat();

  // Histograms.
  cout << nCharged << etaCharged << pTCharged;

  // Done.
  return 0;
}
