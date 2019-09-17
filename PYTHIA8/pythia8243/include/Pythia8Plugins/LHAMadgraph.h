// LHAMadgraph.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Philip Ilten, December 2015.

#ifndef Pythia8_LHAMadgraph_H
#define Pythia8_LHAMadgraph_H

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/JetMatching.h"
#include "Pythia8Plugins/GeneratorInput.h"
#include <unistd.h>
#include <sys/stat.h>

using namespace std;

namespace Pythia8 {

//==========================================================================

// A derived class from LHAup which generates events with MadGraph 5.

// This class automatically generates hard processes with MadGraph 5
// and aMC@NLO, reads in the LHEF file output, and performs matching
// (if requested). For tree-level generation MLM matching is performed
// while FxFx is used for aMC@NLO generation.

// The user can send commands to MadGraph via the readString
// method. Any string begining with "configure " is used for the initial
// MadGraph configuration with "configure " stripped from the
// begining. In general, only the process and run settings need to be
// provided. Run settings must begin with " set" (note the leading
// space). The output and launch commands, random seed, and shower
// choice are automatically handled. For example, the following will
// produce di-muon events from 13 TeV proton proton collisions at NLO
// in QCD:

//    readString("generate p p > mu+ mu- [QCD]")

// The number of events generated per MadGraph run is controlled by
// setEvents, while the random seed is controlled by setSeed. The
// maximum number of jets produced by MadGraph (needed for matching)
// is automatically determined but can be manually specified with
// setJets. In general these methods should not be needed; for further
// details see the method documentation.

// Events are generated with MadGraph utilizing the "gridpack" method
// for MadGraph 5:
//
//    https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/GridDevelopment
//
// and an eqivalent method for aMC@NLO:
//
//    https://answers.launchpad.net/mg5amcnlo/+question/243268
//
// Consequently the run directory ("madgraphrun" by default) does not
// need to be deleted between independent runs with the same
// configuration (excluding random seeds). Indeed, keeping the
// directory significantly speeds the generation process, particularly
// for NLO generation with aMC@NLO as the grid initialization
// can be skipped after the initial run.

class LHAupMadgraph : public LHAup {

public:

  // Types of MadGraph stages.
    enum Stage{Auto, Configure, Generate, Launch};

  // Constructor.
  LHAupMadgraph(Pythia* pythiaIn, bool matchIn = true,
                string dirIn = "madgraphrun", string exeIn = "mg5_aMC");

  // Destructor.
  ~LHAupMadgraph();

  // Read a MadGraph command string.
  bool readString(string line, Stage stage = Auto);

  // Add a MadGraph configuration card to be used.
  void addCard(string src, string dst);

  // Set the number of events to generate per run.
  void setEvents(int eventsIn);

  // Set the random seed and maximum runs.
  bool setSeed(int seedIn, int runsIn = 30081);

  // Set the maximum number of jets produced by MadGraph.
  void setJets(int jetsIn);

  // Set the initialization information.
  bool setInit();

  // Set the event information.
  bool setEvent(int = 0);

protected:

  // Execute a system command.
  bool execute(string line);

  // Write the MadGraph configuration.
  bool configure();

  // Run the generate stage of MadGraph.
  bool generate();

  // Run the launch stage of MadGraph.
  bool launch();

  // Run MadGraph.
  bool run(int eventsIn, int seedIn = -1);

  // Create the LHEF reader.
  bool reader(bool init);

  // The PYTHIA object and LHEF file reader and matching hook.
  Pythia *pythia;
  LHAupLHEF *lhef;
  JetMatchingMadgraph *hook;

  // Stored members.
  int events, seed, runs, nRuns, jets;
  bool match, amcatnlo;
  string dir, exe, lhegz;
  double sigWgt, wgt;
  vector< pair<string, string> > cards;

  // The MadGraph commands for the config, generate, and launch stages.
  vector<string> configureLines, generateLines, launchLines;

  // Vector of whether a command stage has been overridden by the user.
  vector<bool> override;

};

//--------------------------------------------------------------------------

// Constructor.

LHAupMadgraph::LHAupMadgraph(Pythia *pythiaIn, bool matchIn, string dirIn,
                             string exeIn) :
  pythia(pythiaIn), lhef(0), hook(0), events(10000), seed(-1), runs(30081),
  nRuns(0), jets(-1), match(matchIn), dir(dirIn), exe(exeIn),
  lhegz(dirIn + "/events.lhe.gz"), override(vector<bool>(3, false)) {
  mkdir(dir.c_str(), 0777);
  if (pythia) pythia->readString("Beams:frameType = 5");
}

//--------------------------------------------------------------------------

// Destructor.

LHAupMadgraph::~LHAupMadgraph() {if (lhef) delete lhef; if (hook) delete hook;}

//--------------------------------------------------------------------------

// Read a MadGraph command string.

// If the stage is set to Auto, commands beginning with " set" are
// used in the launch stage (these must begin with a single space to
// differentiate from generate stage set commands), commands begining
// with "configure" are used in the configuration stage, and all
// remaining commands (excluding output and launch) are used in the
// generate stage. Output, launch, seed, and shower commands are
// automatically handled. If the user wishes to override commands,
// then the stage can be specified. This will prevent any
// automatically generated commands from being used for that
// stage. This should only be done if the user understands what
// additional commands are needed.

// The system MadGraph configuration will be used if the configuration
// stage is overridden by the user and only blank commands have been
// passed. This is accomplished via
// readString("", LHAupMadgraph::Configure).

bool LHAupMadgraph::readString(string line, Stage stage) {
  if (stage == Auto) {
    if (line.substr(0, 4) == " set") launchLines.push_back(line);
    else if (line.substr(0, 10) == "configure ")
      configureLines.push_back(line.substr(10));
    else if (line.substr(0, 6) != "output" && line.substr(0, 6) != "launch")
      generateLines.push_back(line);
    else return false;
  } else if (stage == Configure) {
    override[Configure] = true; if (line != "") configureLines.push_back(line);
  } else if (stage == Generate) {
    override[Generate] = true; generateLines.push_back(line);
  } else if (stage == Launch) {
    override[Launch] = true; launchLines.push_back(line);
  } else return false;
  return true;
}

//--------------------------------------------------------------------------

// Add a MadGraph configuration card to be used.

// In general, MadGraph should be configured via the readString
// method. However, there are some cases where the user might wish to
// provide an entire configuration card, e.g. setting BSM model
// space. This method allows the user to provide a source card, src,
// which is then copied to <MadGraph run directory>/Cards/dst. These
// cards are copied before any MadGraph processes are launched.

void LHAupMadgraph::addCard(string src, string dst) {
  cards.push_back(make_pair(src, dst));
}

//--------------------------------------------------------------------------

// Set the random seed and maximum allowed runs.

// If the random seed is negative (default of -1), then the MadGraph
// seed is taken as the Pythia parameter "Random:seed", which must be
// greater than 0. If the maximum number of allowed runs is exceeded
// (default of 30081) an error is thrown. The seed for a MadGraph run
// is set as:

//    (random seed - 1) * (maximum runs) + (number of runs) + 1

// MadGraph can only handle random seeds up to 30081 * 30081. So, with
// this strategy, one can generate Pythia jobs with seeds from 1 to
// 30081, with each job running MadGraph less than 30081 times, and
// ensure a fully statistically independent sample. If more than 30081
// jobs are needed, then the maximum allowed runs can be lowered
// accordingly, and if need be, setEvents can be used to increase the
// number of events generated per run.

bool LHAupMadgraph::setSeed(int seedIn, int runsIn) {

  if (!pythia) return false;
  seed = seedIn;
  if (seed < 0) {
    seed = pythia->settings.mode("Random:seed");
    if (seed < 1) {
      pythia->info.errorMsg("Error from LHAupMadgraph::setSeed: the given "
                            "Pythia seed is less than 1."); return false;}
  }
  runs = runsIn;
  if (seed * runs > 30081 * 30081) {
    pythia->info.errorMsg("Error from LHAupMadgraph::setSeed: the given seed "
                          "exceeds the MadGraph limit."); return false;}
  nRuns = 0;
  return true;

}

//--------------------------------------------------------------------------

// Set the number of events to generate per MadGraph run, the default is 10000.

void LHAupMadgraph::setEvents(int eventsIn) {events = eventsIn;}

//--------------------------------------------------------------------------

// Set the number maximum number of jets generated by MadGraph.

// If negative (default of -1) then the number of jets is determined
// automatically. This is the maximum number of jets produced at
// leading order.

void LHAupMadgraph::setJets(int jetsIn) {jets = jetsIn;}

//--------------------------------------------------------------------------

// Execute a system command.

bool LHAupMadgraph::execute(string line) {return system(line.c_str()) != -1;}

//--------------------------------------------------------------------------

// Write the MadGraph configuration.

// If not overridden, the MadGraph configuration is set to prevent the
// output from being opened in a web-browser and stop MadGraph
// updates.

bool LHAupMadgraph::configure() {

  if (override[Configure] && configureLines.size() == 0) return true;
  mkdir((dir + "/.mg5").c_str(), 0777);
  fstream config((dir + "/.mg5/mg5_configuration.txt").c_str(), ios::out);
  for (int iLine = 0; iLine < (int)configureLines.size(); ++iLine)
    config << configureLines[iLine] << "\n";
  if (!override[Configure])
    config << "automatic_html_opening = False\n"
           << "auto_update = 0\n";
  config.close();
  return true;

}

//--------------------------------------------------------------------------

// Run the generate stage of MadGraph.

// The final command of "output <dir> -f -nojpeg\n" is automatically
// set, if not overridden. MadGraph is then run and the output is
// checked. Finally, the configuration is updated and the type of run,
// MadGraph or aMC@NLO, is determined.

bool LHAupMadgraph::generate() {

  // Write the user settings to the generate configuration file.
  if (!pythia) return false;
  fstream config((dir + "/generate.py").c_str(), ios::out);
  for (int iLine = 0; iLine < (int)generateLines.size(); ++iLine)
    config << generateLines[iLine] << "\n";
  if (!override[Generate]) config << "output " << dir << "/tmp -f -nojpeg\n";
  config.close();

  // Run MadGraph and check output.
  fstream orig((dir + "/.mg5/mg5_configuration.txt").c_str(), ios::in);
  char *home = getenv("HOME");
  setenv("HOME", dir.c_str(), 1);
  bool success = execute(exe + " " + dir + "/generate.py");
  setenv("HOME", home, 1);
  if (!success) {orig.close(); return false;}
  else if (access((dir + "/tmp/Cards/run_card.dat").c_str(), F_OK) == -1) {
    pythia->info.errorMsg("Error from LHAupMadgraph::generate: MadGraph "
                          "failed to produce run_card.dat");
    orig.close(); return false;
  } else execute("mv " + dir + "/tmp/* " + dir + "; rmdir " + dir + "/tmp");

  // Update configuration.
  amcatnlo =
    access((dir + "/Cards/amcatnlo_configuration.txt").c_str(), F_OK) != -1;
  if (orig.good()) {
    fstream copy((dir + "/Cards/" + (amcatnlo ? "amcatnlo" : "me5") +
                  "_configuration.txt").c_str(), ios::out);
    copy << orig.rdbuf(); copy.close();
  }
  orig.close();

  // Copy over any user provided configuration cards.
  for (int iCard = 0; iCard < (int)cards.size(); ++iCard) {
    ifstream src((cards[iCard].first).c_str(), ios::binary);
    ofstream dst((dir + "/Cards/" + cards[iCard].second).c_str(), ios::binary);
    dst << src.rdbuf();
  }
  return true;

}

//--------------------------------------------------------------------------

// Run the launch stage of MadGraph.

// The first command "launch ..." is automatically set, and depends on
// whether aMC@NLO is being used.

bool LHAupMadgraph::launch() {

  // Open the launch configuration file and write the settings.
  if (!pythia) return false;
  fstream config((dir + "/launch.py").c_str(), ios::out);
  if (!override[Launch]) {
    config << "launch " << dir << " -n run";
    if (amcatnlo) config << " -p\n" << "set parton_shower PYTHIA8\n"
                         << "set ickkw 3\n" << "set nevents 0\n"
                         << "set req_acc 0.001\n";
    else config << " -s parton\n" << "set ickkw 1\n" << "set gridpack True\n";
  }

  // Write the user settings.
  for (int iLine = 0; iLine < (int)launchLines.size(); ++iLine)
    config << launchLines[iLine] << "\n";
  if (!override[Launch]) config << "done\n";
  config.close();

  // Fix aMC@NLO linking.
  if (amcatnlo) {
    string line = "cd " + dir + "/MCatNLO/lib; LINKS=`ls`; for LINK in $LINKS;"
      " do TARG=`readlink $LINK`; if [[ $TARG = ../* ]]; then "
      "rm $LINK; ln -s ${TARG:3} $LINK; fi; done";
    if (!execute(line)) {
      pythia->info.errorMsg("Error from LHAupMadgraph::launch: failed to "
                            "link aMC@NLO libraries"); return false;}
  }

  // Run MadGraph and create run scripts.
  if (!execute(exe + " " + dir + "/launch.py")) return false;
  if (amcatnlo) {
    if (access((dir + "/SubProcesses/results.dat").c_str(), F_OK) == -1) {
      pythia->info.errorMsg("Error from LHAupMadgraph::launch: aMC@NLO failed "
                            "to produce results.dat"); return false;}
    fstream script((dir + "/run.sh").c_str(), ios::out);
    script << "#!/usr/bin/env bash\n"
           << "sed -i \"s/.*= *nevents/$1 = nevents/g\" ./Cards/run_card.dat\n"
           << "sed -i \"s/.*= *iseed/$2 = iseed/g\" ./Cards/run_card.dat\n"
           << "./bin/generate_events --parton --nocompile --only_generation "
      "--force --name run\n" << "mv Events/run/events.lhe.gz ./\n";
    script.close(); execute("chmod 755 " + dir + "/run.sh");
  } else {
    string gpk = "run_gridpack.tar.gz";
    if (access((dir + "/" + gpk).c_str(), F_OK) == -1) {
      pythia->info.errorMsg("Error from LHAupMadgraph::launch: MadEvent failed"
                            " to produce " + gpk); return false;}
    string line = "cd " + dir + "; tar -xzf " + gpk + "; cd madevent/lib; "
      "LINK=`readlink libLHAPDF.a`; if [[ $LINK = ../* ]]; then "
      "rm libLHAPDF.a; ln -s ../$LINK libLHAPDF.a; fi; cd ../; "
      "./bin/compile dynamic; ./bin/clean4grid";
    if (!execute(line)) {
      pythia->info.errorMsg("Error from LHAupMadgraph::launch: failed to "
                            "compile MadEvent code"); return false;}
  }
  return true;

}

//--------------------------------------------------------------------------

// Run MadGraph.

bool LHAupMadgraph::run(int eventsIn, int seedIn) {

  if (!pythia) return false;
  if (nRuns >= runs) {
    pythia->info.errorMsg("Error from LHAupMadgraph::run: maximum number "
                          "of allowed runs exceeded."); return false;}
  if (access((dir + "/run.sh").c_str(), F_OK) == -1) return false;
  if (seed < 0 && !setSeed(seed, runs)) return false;
  if (seedIn < 0) seedIn = (seed - 1) * runs + nRuns + 1;
  stringstream line;
  line << "cd " + dir + "; ./run.sh " << eventsIn << " " << seedIn;
  if (!amcatnlo) line << "; rm -rf ./madevent/Events/*";
  if (!execute(line.str())) return false;
  if (access(lhegz.c_str(), F_OK) == -1) return false;
  ++nRuns;
  return true;

}

//--------------------------------------------------------------------------

// Create the LHEF reader.

bool LHAupMadgraph::reader(bool init) {

  // Check valid LHE file.
  if (!pythia) return false;
  if (lhef) delete lhef;
  bool setScales(pythia->settings.flag("Beams:setProductionScalesFromLHEF"));
  lhef = new LHAupLHEF(&pythia->info, lhegz.c_str(), NULL, false, setScales);
  if (!lhef->setInit()) {
    pythia->info.errorMsg("Error from LHAupMadgraph::reader: failed to "
                          "initialize the LHEF reader"); return false;}
  if (lhef->sizeProc() != 1) {
    pythia->info.errorMsg("Error from LHAupMadgraph::reader: number of "
                          "processes is not 1"); return false;}

  if (init) {

    // Determine the cross-section (if needed).
    double sig(lhef->xSec(0)), err(lhef->xErr(0));
    if (!amcatnlo) {
      fstream results((dir + "/madevent/SubProcesses/"
                       "run_results.dat").c_str(), ios::in);
      string v; vector<double> vs;
      while (std::getline(results, v, ' ')) vs.push_back(atof(v.c_str()));
      if (vs.size() < 2) {
        pythia->info.errorMsg("Error from LHAupMadgraph::reader: could not "
                              "extract cross-section"); return false;}
      sig = vs[0]; err = vs[1];
    }

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

// Set the initialization information.

// If shower matching has been requested, then the matching is also
// set up depending on the type of MadGraph output.

bool LHAupMadgraph::setInit() {

  // Initialize MadGraph.
  if (!pythia) return false;
  if (access((dir + "/run.sh").c_str(), F_OK) == -1) {
    if (!configure()) {
      pythia->info.errorMsg("Error from LHAupMadgraph::setInit: failed to "
                        "create the MadGraph configuration"); return false;}
    if (!generate()) {
      pythia->info.errorMsg("Error from LHAupMadgraph::setInit: failed to "
                            "generate the MadGraph process"); return false;}
    if (!launch()) {
      pythia->info.errorMsg("Error from LHAupMadgraph::setInit: failed to "
                            "launch the MadGraph process"); return false;}
  } else
    amcatnlo =
      access((dir + "/Cards/amcatnlo_configuration.txt").c_str(), F_OK) != -1;

  // Set up matching if requested.
  if (match) {

    // Load the MadGraph parameters.
    ifstream card((dir + "/Cards/run_card.dat").c_str());
    string str((istreambuf_iterator<char>(card)), istreambuf_iterator<char>());
    MadgraphPar mad;
    mad.parse(str);
    mad.printParams();

    // Extract maximum number of jets.
    int iLine = jets < 0 ? 0 : generateLines.size();
    for (; iLine < (int)generateLines.size(); ++iLine) {
      string line  = generateLines[iLine];
      size_t found = line.find(">");
      if (found == string::npos) continue;
      else line = line.substr(found);
      stringstream sline(line); string p; int n(0);
      while (std::getline(sline, p, ' ') && p != ",")
        if (p == "j" || p == "g" || p == "u" || p == "d" || p == "c" ||
            p == "s" || p == "b" || p == "u~" || p == "d~" || p == "c~" ||
            p == "s~" || p == "b~") ++n;
      if (n > jets) jets = n;
    }

    // Common settings.
    double etaj   = mad.getParam("etaj");
    Settings &set = pythia->settings;
    set.flag("JetMatching:merge", true);
    set.mode("JetMatching:scheme", 1);
    set.flag("JetMatching:setMad", false);
    set.mode("JetMatching:nQmatch", mad.getParamAsInt("maxjetflavor"));
    set.parm("JetMatching:qCut", mad.getParam("ptj"));
    set.parm("JetMatching:etaJetMax", etaj > 0 ? etaj : 100);
    set.mode("JetMatching:nJetMax", jets);
    set.parm("Check:epTolErr", 1e-2);

    // aMC@NLO settings.
    if (amcatnlo) {
      set.parm("JetMatching:coneRadius", mad.getParam("jetradius"));
      set.mode("JetMatching:slowJetPower", mad.getParam("jetalgo"));
      set.parm("JetMatching:qCutME", mad.getParam("ptj"));
      set.mode("JetMatching:jetAlgorithm", 2);
      set.flag("JetMatching:doFxFx", true);
      set.flag("SpaceShower:MEcorrections", false);
      set.parm("TimeShower:pTmaxMatch", 1);
      set.parm("TimeShower:pTmaxFudge", 1);
      set.flag("TimeShower:MEcorrections", false);
      set.flag("TimeShower:globalRecoil", true);
      set.flag("TimeShower:limitPTmaxGlobal", true);
      set.mode("TimeShower:nMaxGlobalRecoil", 1);
      set.mode("TimeShower:globalRecoilMode", 2);

    // MLM tree-level MadGraph settings.
    } else set.parm("JetMatching:clFact", mad.getParam("alpsfact"));

    // Set the matching hook.
    hook = new JetMatchingMadgraph();
    pythia->setUserHooksPtr(hook);
  }

  // Create the LHEF LHAup object and run setInit.
  if (!run(events)) return false;
  if (!reader(true)) return false;
  listInit();
  return true;

}

//--------------------------------------------------------------------------

// Set the event information.

bool LHAupMadgraph::setEvent(int) {

  // Run setEvent from the LHEF object and launch MadGraph if failed.
  if (!pythia) return false;
  if (!lhef) {
    pythia->info.errorMsg("Error from LHAupMadgraph::setEvent: LHAupLHEF "
                          "object not correctly initialized"); return false;}
  if (!lhef->fileFound()) {
    pythia->info.errorMsg("Error from LHAupMadgraph::setEvent: LHEF "
                          "event file was not found"); return false;}
  if (!lhef->setEvent()) {
    if (!run(events)) return false;
    if (!reader(false)) return false;
    lhef->setEvent();
  }

  // Read the event from the LHEF object.
  setProcess(lhef->idProcess(), lhef->weight(), lhef->scale(),
    lhef->alphaQED(), lhef->alphaQCD());
  for (int ip = 1; ip < lhef->sizePart(); ++ip)
    addParticle(lhef->id(ip), lhef->status(ip), lhef->mother1(ip),
                lhef->mother2(ip), lhef->col1(ip), lhef->col2(ip),
                lhef->px(ip), lhef->py(ip), lhef->pz(ip), lhef->e(ip),
                lhef->m(ip), lhef->tau(ip), lhef->spin(ip), lhef->scale(ip));
  setIdX(lhef->id1(), lhef->id2(), lhef->x1(), lhef->x2());
  setPdf(lhef->id1pdf(), lhef->id2pdf(), lhef->x1pdf(), lhef->x2pdf(),
         lhef->scalePDF(), lhef->pdf1(), lhef->pdf2(), lhef->pdfIsSet());
  return true;

}

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_LHAMadgraph_H
