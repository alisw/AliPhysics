// main112.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This test program will generate p-Pb collisions at sqrt(S_NN)=5TeV
// using Angantyr model for Heavy Ion collisions. The analysis will
// divide the event in centrality classes and measure the charged
// pseudo-rapidity distribution in each class in the way described in
// the ATLAS analysis in arXiv:1508.00848 [hep-ex].

#include "Pythia8/Pythia.h"

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;

int main() {

  Pythia pythia;

  // Setup the beams.
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 1000822080"); // The lead ion.
  pythia.readString("Beams:eA = 4000");
  pythia.readString("Beams:eB = 1570");
  pythia.readString("Beams:frameType = 2");

  // Initialize the Angantyr model to fit the total and semi-includive
  // cross sections in Pythia within some tolerance.
  pythia.readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
  // These parameters are typicall suitable for sqrt(S_NN)=5TeV
  pythia.readString("HeavyIon:SigFitDefPar = "
                    "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
  // A simple genetic algorithm is run for 20 generations to fit the
  // parameters.
  pythia.readString("HeavyIon:SigFitNGen = 20");

  // There will be eight centrality bins based on the sum transverse
  // emergy in a rapidity interval between -4.9 and -3.2. The borders
  // between the classes have been read off the plot in the paper:
  double explim[] = {90.0, 66.0, 53.0, 41.0, 32.0, 24.0, 13.0, 6.0};
  // Alternatively we can obtain the borders from the generated
  // transverse energy spectrum. The default settings should give
  // approximately the following:
  double genlim[] = {77.5, 54.5, 44.4, 34.1, 27.6, 22.5, 14.2, 3.5};
  // If you change any parameters these should also be changed.

  // The upper edge of the correponding percentiles:
  double pclim[] = {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.9};

  // Book the pseudorapidity histograms and get counter for sum of
  // event weights:
  typedef map<double,int,std::greater<double> > MapIdx;
  MapIdx expetaidx, genetaidx;
  vector<Hist*> expetadist(8), genetadist(8);
  string expetaname("EtadistCexp"), genetaname("EtadistCgen");
  vector<double> expsumw(8, 0.0), gensumw(8, 0.0);
  for ( int i = 0; i < 8; ++i ) {
    expetaidx[explim[i]] = i;
    expetadist[i] = new Hist(expetaname + char('0' + i), 54, -2.7, 2.7);
    genetaidx[genlim[i]] = i;
    genetadist[i] = new Hist(genetaname + char('0' + i), 54, -2.7, 2.7);
  }

  // Book histogram for the centrality measure.
  Hist sumet("SumETfwd", 100, 0.0, 200.0);

  // Also make a map of all weight to check the generated centrality
  // classes.
  multimap<double,double> gencent;

  // Book a histogram for the distribution of number of wounded
  // nucleons.
  Hist wounded("Nwounded", 60, -0.5, 59.5);

  // Sum up the weights of all generated events.
  double sumw = 0.0;

  // Initialise Pythia.
  pythia.init();

  // Loop over events.
  int nEvents = 1000;
  for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {
    if ( !pythia.next() ) continue;

    // First sum up transverse energy for centrality measure and also
    // check that the trigger requiring ar least one charged particle
    // forward and backward.
    double etfwd = 0.0;
    bool trigfwd = false;
    bool trigbwd = false;
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() ) {
        double eta = p.eta();
        if ( p.isCharged() && p.pT() > 0.1 && eta < -2.09 && eta > -3.84 )
          trigfwd = true;
        if ( p.isCharged() && p.pT() > 0.1 && eta > 2.09 && eta < 3.84 )
          trigbwd = true;
        if ( p.pT() > 0.1 && eta < -3.2 && eta > -4.9 )
          etfwd += p.eT();
      }
    }
    // Skip if not triggered
    if ( !(trigfwd && trigbwd) ) continue;

    // Keep track of the sum of waights
    double weight = pythia.info.weight();
    sumw += weight;

    // Histogram and save the summed Et.
    sumet.fill(etfwd, weight);
    gencent.insert(make_pair(etfwd, weight));

    // Also fill the number of (absorptively and diffractively)
    // wounded nucleaons.
    int nw = pythia.info.hiinfo->nAbsTarg() +
             pythia.info.hiinfo->nDiffTarg();
    wounded.fill(nw, weight);

    // Find the correct centrality histograms.
    MapIdx::iterator expit =  expetaidx.upper_bound(etfwd);
    int expidx = expit== expetaidx.end()? -1: expit->second;
    MapIdx::iterator genit = genetaidx.upper_bound(etfwd);
    int genidx = genit== genetaidx.end()? -1: genit->second;

    // Sum the weights in the centrality classes, skip if not in a class.
    if ( expidx < 0 && genidx < 0 ) continue;
    if ( expidx >= 0 ) expsumw[expidx] += weight;
    if ( genidx >= 0 ) gensumw[genidx] += weight;

    // Go through the event again and fill the eta distributions.
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() && p.isCharged() &&
           abs(p.eta()) < 2.7 && p.pT() > 0.1 ) {
        if ( expidx >= 0 ) expetadist[expidx]->fill(p.eta(), weight);
        if ( genidx >= 0 ) genetadist[genidx]->fill(p.eta(), weight);
      }
    }
  }

  // The run is over, so we write out some statistics.


  // Now, we just have to normalize and prtint out the histograms. We
  // choose to print the histograms to a file that can be read by
  // eg. gnuplot.
  ofstream ofs("main112.dat");

  sumet /= sumw*2.0;
  ofs << "# " << sumet.getTitle() << endl;
  sumet.table(ofs);

  wounded /= sumw;
  ofs << "\n# " << wounded.getTitle() << endl;
  wounded.table(ofs);


  // Print out the centrality binned eta distributions and delete the
  // heap-allocate histograms.
  for ( int idx = 0; idx < 8; ++idx ) {
    *expetadist[idx] /= expsumw[idx]*0.1;
    ofs << "\n# " << expetadist[idx]->getTitle() << endl;
    expetadist[idx]->table(ofs);
    delete expetadist[idx];
    *genetadist[idx] /= gensumw[idx]*0.1;
    ofs << "\n# " << genetadist[idx]->getTitle() << endl;
    genetadist[idx]->table(ofs);
    delete genetadist[idx];
  }

  // Befor we end, we want to check that our generated centrality
  // classes were the same as we guessed.
  double curr = 0.0;
  double prev = 0.0;
  double acc = 0.0;
  int idxa = 7;
  double lim = sumw*(1.0 - pclim[idxa]);
  vector<double> newlim(8);
  for ( multimap<double, double>::iterator it = gencent.begin();
        it != gencent.end(); ++it ) {
    prev = curr;
    curr = it->first;
    double w = it->second;
    if ( acc < lim && acc + w >= lim ) {
      newlim[idxa--] = prev + (curr - prev)*(lim - acc)/w;
      lim = sumw*(1.0 - pclim[idxa]);
    }
    acc += w;
  }

  cout << "The generated limits between centrality classes in this run:\n"
       << "   %   assumed    actual      data\n";
  for ( int idx = 0; idx < 8; ++idx )
    cout << setw(4) << int(pclim[idx]*100.0 + 0.5)
         << setw(10) << fixed << setprecision(1) << genlim[idx]
         << setw(10) << fixed << setprecision(1) << newlim[idx]
         << setw(10) << fixed << setprecision(1) << explim[idx] << endl;

  // And we're done!
  return 0;
}
