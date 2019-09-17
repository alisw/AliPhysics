// Pythia8Rivet.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#ifndef PYTHIA8RIVET_H
#define PYTHIA8RIVET_H

#include "Pythia8/HIUserHooks.h"
#include "Pythia8Plugins/HepMC2.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include <set>
#include <vector>
#include <string>
#include "Rivet/Rivet.hh"

namespace Pythia8 {

using namespace std;

/**
 * Simplified interface to the Rivet program. Remember to link with
 * pythia and -lhepmcinterface -lHepMC -lRivet
 *
 * Usage: (1) Create an object giving the pythia object and a filename
 * as arguments. (2) Repeatedly specify (the name of an) analysis with
 * the addAnalysis() function, possibly with analysis parameters.
 * (3) initialize the underlying Rivet  object with the init() function.
 * (4) Analyze an event with the operator() function. (5) Dump the
 * histograms to a file with the done() function.
 */
class Pythia8Rivet {

public:

  /**
   * The constructor needs to have the main Pythia object, and the
   * name of the file where the histograms are dumped.
   */
  Pythia8Rivet(Pythia & pytin, string fname)
    : pythia(&pytin), filename(fname), rivet(0), igBeam(false) {}

  /**
   * The destructor will write the histogram file if this has not
   * already been done.
   */
  ~Pythia8Rivet() {
    done();
  }

  /**
   * Add the name of an analysis to be performed, with a list
   * of analysis parameters.
   */
  void addAnalysis(string ana) {
    analyses.insert(ana);
  }

  /**
   * Add a YODA file pre-load pre-filled histogram from.
   */
  void addPreload(string prel) {
    preloads.push_back(prel);
  }

  /**
   * Set "ignore beams" flag.
   */
  void ignoreBeams(bool flagIn) {
    igBeam = flagIn;
  }

  /**
   * Add an (optional) run name for Rivet internal use.
   */
  void addRunName(const string runname) {
    rname = runname;
  }
  /**
   * Initialize Rivet. Will do nothing if Rivet was already initialized
   */
  void init(const HepMC::GenEvent & gev) {
    if ( rivet ) return;
    rivet = new Rivet::AnalysisHandler(rname);
    rivet->setIgnoreBeams(igBeam);
    Rivet::addAnalysisLibPath(".");
    for(int i = 0, N = preloads.size(); i < N; ++i)
      rivet->readData(preloads[i]);
    for (set<string>::iterator it = analyses.begin();
      it != analyses.end(); ++it) {
      rivet->addAnalysis(*it);
    }
    rivet->init(gev);
  }

  /**
   * Analyze the default Pythia event. Will do nothing if Rivet has
   * not been intialized.
   */
  void operator()() {
    this->operator()(pythia->event, -1, &pythia->info, &pythia->settings);
  }

  /**
   * Analyze the given event.
   */
  void operator()(Event & event, int ievnum = -1, Pythia8::Info* pyinfo = 0,
                  Pythia8::Settings* pyset = 0, bool append = false,
                  HepMC::GenParticle* rootParticle = 0, int iBarcode = -1) {
    HepMC::GenEvent geneve;
    converter.fill_next_event(event, &geneve, ievnum, pyinfo, pyset,
                              append, rootParticle, iBarcode);
    if ( pyinfo && pyinfo->hiinfo ) {
      HepMC::HeavyIon ion;
      ion.set_Ncoll_hard(pyinfo->hiinfo->nCollNDTot());
      ion.set_Ncoll(pyinfo->hiinfo->nAbsProj() +
                    pyinfo->hiinfo->nDiffProj() +
                    pyinfo->hiinfo->nAbsTarg() +
                    pyinfo->hiinfo->nDiffTarg() -
                    pyinfo->hiinfo->nCollND() -
                    pyinfo->hiinfo->nCollDD());
      ion.set_Npart_proj(pyinfo->hiinfo->nAbsProj() +
                         pyinfo->hiinfo->nDiffProj());
      ion.set_Npart_targ(pyinfo->hiinfo->nAbsTarg() +
                         pyinfo->hiinfo->nDiffTarg());
      ion.set_impact_parameter(pyinfo->hiinfo->b());
      geneve.set_heavy_ion(ion);
    }
    if ( !rivet ) init(geneve);
    rivet->analyze(geneve);
  }

  /**
   * Writes histograms to file and deletes the Rivet object. Does
   * nothing if Rivet was not initialized.
   */
  void done() {
    if ( !rivet ) return;
    rivet->finalize();
    rivet->writeData(filename);
    delete rivet;
    rivet = 0;
  }

private:

  /**
   * The main pythia object.
   */
  Pythia * pythia;

  /**
   * The name of the file where the histograms are dumped.
   */
  string filename;

  /**
   * Analyses with optional analysis parameters.
   */
  set<string> analyses;

  /**
   * The names of YODA files to preload.
   */
  vector<string> preloads;

  /**
   * The Rivet object.
   */
  Rivet::AnalysisHandler * rivet;

  /**
   * The HepMC converter
   */
  HepMC::Pythia8ToHepMC converter;

  /**
   * The Rivet run name
   */
  string rname;

  /**
   * Ignore beams flag.
   */
  bool igBeam;

};

}

#endif
