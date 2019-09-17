// main113.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This test program will generate Pb-Pb collisions at
// sqrt(S_NN)=2.76TeV using the Angantyr model for Heavy Ion
// collisions. The analysis will divide the event in centrality
// classes using the same observable as was used for p-Pb in the ATLAS
// analysis in arXiv:1508.00848 [hep-ex] (see main112.cc). The
// centrality classes are same as in the ALICE analysis in
// arXiv:1012.1657 [nucl-ex] although the actual observable used is
// not the same. Histograms of multiplicity distributions are measured
// for each centrality percentile.

// Note that heavy ion collisions are computationally quite CPU
// intensive and generating a single event will take around a second
// on a reasonable desktop. To get reasonable statistics, this program
// will take a couple of hours to run.

#include "Pythia8/Pythia.h"

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;

int main() {

  Pythia pythia;

  // Setup the beams.
  pythia.readString("Beams:idA = 1000822080");
  pythia.readString("Beams:idB = 1000822080"); // The lead ion.
  pythia.readString("Beams:eCM = 2760.0");
  pythia.readString("Beams:frameType = 1");

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

  // There will be nine centrality bins based on the sum transverse
  // emergy in a rapidity interval between 3.2 and 4.9 obtained from
  // the borders from the generated transverse energy spectrum. The
  // default settings should give approximately the following:
  double genlim[] = {2979.4, 2400.1, 1587.5, 1028.8, 669.9,
                     397.4, 220.3, 116.3, 54.5};
  // If you change any parameters these should also be changed.

  // The upper edge of the correponding percentiles:
  double pclim[] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  // Book the pseudorapidity and multiplicity histograms and get
  // counters for number of events and sum of event weights:
  typedef map<double,int,std::greater<double> > MapIdx;
  MapIdx genetaidx;
  vector<Hist*> etadist(9), lmult(9), hmult(9);
  string etaname("EtadistC"), mname("MultC");
  vector<double> gensumw(9, 0.0), gensumn(9, 0.0);
  for ( int i = 0; i < 9; ++i ) {
    genetaidx[genlim[i]] = i;
    etadist[i] = new Hist(etaname + char('0' + i), 54, -2.7, 2.7);
    hmult[i] = new Hist(mname + 'H' + char('0' + i),
                        75, -0.5, 2999.5);
    lmult[i] = new Hist(mname + 'L' + char('0' + i),
                        75, -0.5, 299.5);
  }

  // Book histogram for the centrality measure.
  Hist sumet("SumETfwd", 200, 0.0, 4000.0);

  // Also make a map of all weight to check the generated centrality
  // classes.
  multimap<double,double> gencent;

  // Book a histogram for the distribution of number of wounded
  // nucleons.
  Hist wounded("Nwounded", 209, -0.5, 417.5);

  // Profile for average central multiplicity and number of wounded
  // nucleons as a function of centrality (with errors).
  vector<double> cmult(9, 0.0), cmult2(9, 0.0);
  vector<double> wound(9, 0.0), wound2(9, 0.0);

  // Sum up the weights of all generated events.
  double sumw = 0.0;

  // Initialise Pythia.
  pythia.init();

  // Loop over events.
  int nEvents = 10000;
  for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {
    if ( !pythia.next() ) continue;

    // First sum up transverse energy for centrality measure and also
    // check that the trigger requiring ar least one charged particle
    // forward and backward.
    double etfwd = 0.0;
    bool trigfwd = false;
    bool trigbwd = false;
    int nc = 0;
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() ) {
        double eta = p.eta();
        if ( p.isCharged() && p.pT() > 0.1 && eta < -2.09 && eta > -3.84 )
          trigfwd = true;
        if ( p.isCharged() && p.pT() > 0.1 && eta > 2.09 && eta < 3.84 )
          trigbwd = true;
        if ( p.pT() > 0.1 && abs(eta) > 3.2 && abs(eta) < 4.9 )
          etfwd += p.eT();
        if ( p.isCharged() && p.pT() > 0.1 && abs(eta) < 0.5 ) ++nc;
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
             pythia.info.hiinfo->nDiffTarg() +
             pythia.info.hiinfo->nAbsProj() +
             pythia.info.hiinfo->nDiffProj();
    wounded.fill(nw, weight);

    // Find the correct centrality histograms.
    MapIdx::iterator genit = genetaidx.upper_bound(etfwd);
    int genidx = genit== genetaidx.end()? -1: genit->second;

    // Sum the weights in the centrality classes, skip if not in a class.
    if ( genidx < 0 ) continue;
    gensumw[genidx] += weight;
    hmult[genidx]->fill(nc, weight);
    lmult[genidx]->fill(nc, weight);
    gensumn[genidx] += 1.0;
    cmult[genidx] += nc*weight;
    cmult2[genidx] += nc*nc*weight;
    wound[genidx] += nw*weight;
    wound2[genidx] += nw*nw*weight;

    // Go through the event again and fill the eta distributions.
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() && p.isCharged() &&
           abs(p.eta()) < 2.7 && p.pT() > 0.1 ) {
         etadist[genidx]->fill(p.eta(), weight);
      }
    }
  }

  // The run is over, so we write out some statistics.


  // Now, we just have to normalize and prtint out the histograms. We
  // choose to print the histograms to a file that can be read by
  // eg. gnuplot.
  ofstream ofs("main113.dat");

  sumet /= sumw*2.0;
  ofs << "# " << sumet.getTitle() << endl;
  sumet.table(ofs);

  wounded /= sumw*2.0;
  ofs << "\n# " << wounded.getTitle() << endl;
  wounded.table(ofs);

  // Print out the centrality binned eta distributions and delete the
  // heap-allocate histograms.
  for ( int idx = 0; idx < 9; ++idx ) {
    *hmult[idx] /= gensumw[idx]*40.0;
    ofs << "\n# " << hmult[idx]->getTitle() << endl;
    hmult[idx]->table(ofs);
    delete hmult[idx];
    *lmult[idx] /= gensumw[idx]*4.0;
    ofs << "\n# " << lmult[idx]->getTitle() << endl;
    lmult[idx]->table(ofs);
    delete lmult[idx];
    *etadist[idx] /= gensumw[idx]*0.1;
    ofs << "\n# " << etadist[idx]->getTitle() << endl;
    etadist[idx]->table(ofs);
    delete etadist[idx];
  }

  // Print out average central charged multiplicity as a function of
  // centrality.
  ofs << "\n# Nch0\n";
  for ( int idx = 0; idx < 9; ++idx ) {
    double Nch = cmult[idx]/gensumw[idx];
    cmult2[idx] = (cmult2[idx]/gensumw[idx] - pow2(Nch))/gensumn[idx];
    ofs << setprecision(2) << setw(4) << int(pclim[idx]*100.0 + 0.5)
        << setw(10) << Nch << setw(10) << sqrt(cmult2[idx]) <<endl;
  }
  ofs << "\n# Nwc\n";
  for ( int idx = 0; idx < 9; ++idx ) {
    double Nw = wound[idx]/gensumw[idx];
    wound2[idx] = (wound2[idx]/gensumw[idx] - pow2(Nw))/gensumn[idx];
    ofs << setprecision(2) << setw(4) << int(pclim[idx]*100.0 + 0.5)
        << setw(10) << Nw << setw(10) << sqrt(wound2[idx]) <<endl;
  }

  // Befor we end we print out some statistics. Also, we want to check
  // that our generated centrality classes were the same as we
  // guessed.
  pythia.stat();
  double curr = 0.0;
  double prev = 0.0;
  double acc = 0.0;
  int idxa = 8;
  double lim = sumw*(1.0 - pclim[idxa]);
  vector<double> newlim(9);
  for ( multimap<double, double>::iterator it = gencent.begin();
        it != gencent.end(); ++it ) {
    prev = curr;
    curr = it->first;
    double w = it->second;
    if ( acc < lim && acc + w >= lim ) {
      newlim[idxa--] = prev + (curr - prev)*(lim - acc)/w;
      if ( idxa < 0 ) break;
      lim = sumw*(1.0 - pclim[idxa]);
    }
    acc += w;
  }

  cout << "The generated limits between centrality classes in this run:\n"
       << "   %   assumed    actual\n";
  for ( int idx = 0; idx < 9; ++idx )
    cout << setw(4) << int(pclim[idx]*100.0 + 0.5)
         << setw(10) << fixed << setprecision(1) << genlim[idx]
         << setw(10) << fixed << setprecision(1) << newlim[idx] << endl;

  // And we're done!
  return 0;
}
