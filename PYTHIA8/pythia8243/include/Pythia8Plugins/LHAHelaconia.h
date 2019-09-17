// LHAHelaconia.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Philip Ilten, December 2017.

#ifndef Pythia8_LHAHelaconia_H
#define Pythia8_LHAHelaconia_H

#include "Pythia8/Pythia.h"
#include <unistd.h>
#include <sys/stat.h>
#include <limits>

namespace Pythia8 {

//==========================================================================

// A derived class from LHAup which generates events with HelacOnia.

// This class automatically generates hard processes with HelacOnia
// and reads in the LHEF file output.

// The user can send commands to HelacOnia via the readString
// method. The launch command, random seed, and shower choice are
// automatically handled. For example, the following will produce
// J/psi events from 13 TeV proton proton collisions:

//    readString("generate u u~ > cc~(3S11) g")

// The number of events generated per HelacOnia run is controlled by
// setEvents, while the random seed is controlled by setSeed.

class LHAupHelaconia : public LHAup {

public:

  // Constructor.
  LHAupHelaconia(Pythia *pythiaIn, string dirIn = "helaconiarun",
                 string exeIn = "ho_cluster");

  // Destructor.
  ~LHAupHelaconia();

  // Read a HelacOnia command string.
  bool readString(string line);

  // Set the number of events to generate per run.
  void setEvents(int eventsIn);

  // Set the random seed and maximum runs.
  bool setSeed(int seedIn, int runsIn = 30081);

  // Set the initialization information.
  bool setInit();

  // Set the event information.
  bool setEvent(int = 0);

protected:

  // Execute a system command.
  bool execute(string line);

  // Run HelacOnia.
  bool run(int eventsIn, int seedIn = -1);

  // Create the LHEF reader.
  bool reader(bool init);

  // Convert the color octet HelacOnia ID to a Pythia 8 ID.
  int convert(int idIn);

  // The PYTHIA object and LHEF file reader and matching hook.
  Pythia *pythia;
  LHAupLHEF *lhef;

  // Stored members.
  int events, seed, runs, nRuns, nId, nQ, nR, nL, nJ;
  string dir, exe, lhegz;
  double sigWgt, wgt, mQ;

  // The HelacOnia commands.
  vector<string> lines;

};

//--------------------------------------------------------------------------

// Constructor.

LHAupHelaconia::LHAupHelaconia(Pythia *pythiaIn, string dirIn, string exeIn) :
  pythia(pythiaIn), lhef(0), events(10000), seed(-1), runs(30081),
  nRuns(0), nId(443), nQ(4), nR(0), nL(0), nJ(3),
  dir(dirIn), exe(exeIn), lhegz(dirIn + "/events.lhe"), mQ(-2) {
  mkdir(dir.c_str(), 0777);
  if (pythia) pythia->readString("Beams:frameType = 5");
  pythia->settings.addMode("Onia:state", -1, false, false, 0, 0);
}

//--------------------------------------------------------------------------

// Destructor.

LHAupHelaconia::~LHAupHelaconia() {if (lhef) delete lhef;}

//--------------------------------------------------------------------------

// Read a HelacOnia command string.

// The special command "set state = <pid>" is not passed to HelacOnia,
// but rather is used to set the color singlet state being produced
// from the color octet state (if a color octet state is being
// produced). If color octet production is enabled then the
// appropriate quark mass is modified to half the mass of the color
// octet state plus the color octet mass splitting.

bool LHAupHelaconia::readString(string line) {

  size_t n = line.find("state");
  if (line.find("8)") != string::npos) mQ = -1;
  if (n != string::npos && pythia) {
    pythia->settings.readString("Onia:" + line.substr(n));
    nId = abs(pythia->settings.mode("Onia:state"));
    nQ  = int(nId/1e2) % 10;
    nR  = int(nId/1e5) % 10;
    nL  = int(nId/1e4) % 10;
    nJ  = int(nId/1e0) % 10;
  } else lines.push_back(line);
  return true;

}

//--------------------------------------------------------------------------

// Set the random seed and maximum allowed runs.

// If the random seed is negative (default of -1), then the HelacOnia
// seed is taken as the Pythia parameter "Random:seed", which must be
// greater than 0. If the maximum number of allowed runs is exceeded
// (default of 30081) an error is thrown. The seed for a HelacOnia run
// is set as:

//    (random seed - 1) * (maximum runs) + (number of runs) + 1

// HelacOnia can only handle random seeds up to 30081 * 30081. So, with
// this strategy, one can generate Pythia jobs with seeds from 0 to
// 30081, with each job running HelacOnia less than 30081 times, and
// ensure a fully statistically independent sample. If more than 30081
// jobs are needed, then the maximum allowed runs can be lowered
// accordingly, and if need be, setEvents can be used to increase the
// number of events generated per run.

bool LHAupHelaconia::setSeed(int seedIn, int runsIn) {

  if (!pythia) return false;
  seed = seedIn;
  if (seed < 0) {
    seed = pythia->settings.mode("Random:seed");
    if (seed < 1) {
      pythia->info.errorMsg("Error from LHAupHelaconia::setSeed: the given "
                            "Pythia seed is less than 1."); return false;}
  }
  runs = runsIn;
  if (seed * runs > 30081 * 30081) {
    pythia->info.errorMsg("Error from LHAupHelaconia::setSeed: the given seed "
                          "exceeds the HelacOnia limit."); return false;}
  nRuns = 0;
  return true;

}

//--------------------------------------------------------------------------

// Set the number of events to generate per HelacOnia run; default is 10000.

void LHAupHelaconia::setEvents(int eventsIn) {events = eventsIn;}

//--------------------------------------------------------------------------

// Execute a system command.

bool LHAupHelaconia::execute(string line) {return system(line.c_str()) != -1;}

//--------------------------------------------------------------------------

// Run HelacOnia.

bool LHAupHelaconia::run(int eventsIn, int seedIn) {

  // Set up run and seed.
  if (!pythia) return false;
  if (nRuns >= runs) {
    pythia->info.errorMsg("Error from LHAupHelaconia::run: maximum number "
                          "of allowed runs exceeded."); return false;}
  if (seed < 0 && !setSeed(seed, runs)) return false;
  if (seedIn < 0) seedIn = (seed - 1) * runs + nRuns + 1;

  // Determine the heavy quark mass.
  if (mQ == -1)
    mQ = (pythia->particleData.m0(nId)
          + pythia->settings.parm("Onia:massSplit"))/2.0;

  // Write the generation file.
  if (!pythia) return false;
  fstream config((dir + "/generate.py").c_str(), ios::out);
  for (int iLine = 0; iLine < (int)lines.size(); ++iLine)
    config << lines[iLine] << "\n";
  config << "set seed = " << seedIn << "\n"
         << "set unwgt = T\n"
         << "set unwevt = " << eventsIn << "\n"
         << "set preunw = " << 3.0/2.0*eventsIn << "\n";
  if (mQ > 0) config << "set " << (nQ == 4 ? "c" : "b") << "mass = " << mQ
                     << "\n";
  config << "launch\n";
  config.close();

  // Create the event shuffler.
  fstream shuffle((dir + "/shuffle.py").c_str(), ios::out);
  shuffle <<
    "import random, os\n"
    "random.seed(" << seedIn << "); tag, pre, post, events = '', [], [], []\n"
    "for line in open('events.lhe').readlines():\n"
    "    if line.strip().startswith('<'):\n"
    "        tag = line.strip()\n"
    "        if tag == '<event>':  events += ['<event>\\n']; continue\n"
    "        if tag == '</event>': events[-1] += '</event>\\n'; continue\n"
    "    if tag == '<event>': events[-1] += line\n"
    "    elif len(events) == 0: pre += [line]\n"
    "    else: post += [line]\n"
    "random.shuffle(events); os.unlink('events.lhe')\n"
    "open('events.lhe', 'w').writelines(pre + events + post)\n";
  shuffle.close();

  // Execute the run commands.
  if (!execute("rm -rf " + dir + "/PROC* " + lhegz)) return false;
  if (!execute("cd " + dir + "; cat generate.py | " + exe)) return false;
  if (!execute("cd " + dir + "; ln -s PROC_HO_0/P0_calc_0/output/*.lhe "
               "events.lhe;# python shuffle.py")) return false;
  if (access(lhegz.c_str(), F_OK) == -1) return false;
  ++nRuns;
  return true;

}

//--------------------------------------------------------------------------

// Create the LHEF reader.

bool LHAupHelaconia::reader(bool init) {

  // Check valid LHE file.
  if (!pythia) return false;
  if (lhef) delete lhef;
  bool setScales(pythia->settings.flag("Beams:setProductionScalesFromLHEF"));
  lhef = new LHAupLHEF(&pythia->info, lhegz.c_str(), NULL, false, setScales);
  if (!lhef->setInit()) {
    pythia->info.errorMsg("Error from LHAupHelaconia::reader: failed to "
                          "initialize the LHEF reader"); return false;}
  if (lhef->sizeProc() != 1) {
    pythia->info.errorMsg("Error from LHAupHelaconia::reader: number of "
                          "processes is not 1"); return false;}

  if (init) {

    // Determine the cross-section (if needed).
    double sig(lhef->xSec(0)), err(lhef->xErr(0));

    // Set the info.
    setBeamA(lhef->idBeamA(), lhef->eBeamA(), lhef->pdfGroupBeamA(),
             lhef->pdfSetBeamA());
    setBeamB(lhef->idBeamB(), lhef->eBeamB(), lhef->pdfGroupBeamB(),
             lhef->pdfSetBeamB());
    setStrategy(lhef->strategy());
    addProcess(lhef->idProcess(0), sig, err, lhef->xMax(0));
    xSecSumSave = sig; xErrSumSave = err;
  }
  return true;

}

//--------------------------------------------------------------------------

// Convert the color octet HelacOnia ID to a Pythia 8 ID.

int LHAupHelaconia::convert(int idIn) {

    if (abs(idIn) < 9900000) return idIn;
    else idIn = abs(idIn) - 9900000;
    int nS = 2;
    if (idIn == nQ*110 + 3) nS = 0;
    else if (idIn == nQ*110 + 1) nS = 1;
    return 9900000 + 10000*nQ + 1000*nS + 100*nR + 10*nL + nJ;

}

//--------------------------------------------------------------------------

// Set the initialization information.

bool LHAupHelaconia::setInit() {

  // Create the LHEF LHAup object and run setInit.
  if (!pythia) return false;
  if (!run(events)) return false;
  if (!reader(true)) return false;
  listInit();
  return true;

}

//--------------------------------------------------------------------------

// Set the event information.

bool LHAupHelaconia::setEvent(int) {

  // Run setEvent from the LHEF object and launch HelacOnia if failed.
  if (!pythia) return false;
  if (!lhef) {
    pythia->info.errorMsg("Error from LHAupHelaconia::setEvent: LHAupLHEF "
                          "object not correctly initialized"); return false;}
  if (!lhef->fileFound()) {
    pythia->info.errorMsg("Error from LHAupHelaconia::setEvent: LHEF "
                          "event file was not found"); return false;}
  if (!lhef->setEvent()) {
    if (!run(events)) return false;
    if (!reader(false)) return false;
    lhef->setEvent();
  }

  // Read the event from the LHEF object.
  particlesSave.clear();
  int mom1, mom2;
  for (int ip = 1; ip < lhef->sizePart(); ++ip) {
    mom1 = lhef->mother1(ip);
    mom2 = lhef->mother2(ip);
    particlesSave.push_back
      (LHAParticle(convert(lhef->id(ip)),
                   lhef->status(ip), mom1, mom2, lhef->col1(ip),
                   lhef->col2(ip), lhef->px(ip), lhef->py(ip), lhef->pz(ip),
                   lhef->e(ip), lhef->m(ip), lhef->tau(ip), lhef->spin(ip),
                   lhef->scale(ip)));
    if (mom1 > 0 && mom1 < (int)particlesSave.size() && mom2 == 0)
      particlesSave[mom1 - 1].statusPart = 2;
  }

  // Write the event.
  setProcess(lhef->idProcess(), lhef->weight(), lhef->scale(),
    lhef->alphaQED(), lhef->alphaQCD());
  for (int ip = 0; ip < (int)particlesSave.size(); ++ip)
    addParticle(particlesSave[ip]);
  setIdX(lhef->id1(), lhef->id2(), lhef->x1(), lhef->x2());
  setPdf(lhef->id1pdf(), lhef->id2pdf(), lhef->x1pdf(), lhef->x2pdf(),
         lhef->scalePDF(), lhef->pdf1(), lhef->pdf2(), lhef->pdfIsSet());
  return true;

}

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_LHAHelaconia_H
