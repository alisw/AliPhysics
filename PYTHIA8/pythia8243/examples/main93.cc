// main93.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"

#include "Pythia8Plugins/HepMC2.h"
#include "Pythia8Plugins/Pythia8Rivet.h"
#ifdef USE_ROOT
#include "TTree.h"
#include "TFile.h"
#include "main93.h"
#endif

using namespace Pythia8;

// Helper class to parse command line options.
class InputParser {

public:

  InputParser (int &argc, char **argv) {
    for (int i = 1; i < argc; ++i) arglist.push_back(string(argv[i]));
  }

  const string& getOption(const string &opt) const {
    vector<string>::const_iterator itr = find(arglist.begin(),
      arglist.end(), opt);
    if (itr != arglist.end() && ++itr != arglist.end()) return *itr;
    return "";
  }

  bool hasOption(const string &opt) const {
    return find(arglist.begin(), arglist.end(), opt) != arglist.end();
  }

private:
  vector<string> arglist;
};


int main(int argc, char* argv[]) {

  // Parser object for command line input.
  InputParser ip(argc, argv);

  // Print help text and exit.
  if(ip.hasOption("-h") || ip.hasOption("--help")) {
    cout << "Usage: Run Pythia with cmnd file input, and get Rivet, HepMC or\n"
      "standard Pythia output.\n" << endl;
    cout << "Examples:\n\n\t ./main93 [options] \n\n\t or\n\n\t ./main93 -c "
      "main93.cmnd -n 1000 -o myoutput\n" << endl;
    cout << "Options:\n"
      "\t -h, --help\n\t\t Show this help message and exit.\n"
      "\t -c CMND-FILE\n\t\t Use this user-written command file.\n"
      "\t -c2 CMND-FILE2\n\t\t Use a second cmnd file, loaded after "
       "the first.\n"
      "\t \t Useful for eg. tuning studies.\n"
      "\t -s SEED \n\t\t Specify seed for the random number generator.\n"
      "\t -o OUT \n\t\t Specify output prefix. Rivet histograms becomes "
      "OUTRivet.yoda.\n"
      "\t -n NEVENTS\n\t\t Number of events.\n"
      "\t -l \n\t\t Silence the splash screen.\n"
        << endl;
     cout << "Additional options in cmnd file:\n"
       "A few extra commands can be added to the cmnd file, compared "
       "to normal.\n"
       "\t Main:runRivet = on \n\t\tRun Rivet analyses (requires a\n"
       "\t\tworking installation of Rivet, linked to main93).\n"
       "\t Main:analyses = ANALYSIS1,ANALYSIS2,...\n "
       "\t\tA comma separated list of desired Rivet analyses to be run.\n"
       "\t\tAnalyses can be post-fixed with Rivet analysis parameters:\n"
        "\t\tANALYSIS:parm->value.\n"
       "\t Main:rivetRunName = STRING \n\t\tAdd an optional run name to\n"
       "\t\tthe Rivet analysis.\n"
       "\t Main:rivetIgnoreBeams = on\n\t\t Ignore beams in Rivet. \n"
       "\t Main:writeHepMC = on \n\t\tWrite HepMC output (requires\n"
       "\t\ta working installation of HepMC, linked to main93.)\n"
       "\t Main:writeRoot = on \n\t\tWrite a root tree defined in the\n"
       "\t\tmain93.h header file. Requires a working installation of Root, \n"
       "\t\tproperly linked to Pythia.\n"
       "\t Main:outputLog = on\n\t\tRedirect output to a logfile called\n"
       "\t\tpythia.log. Can be prefixed with -o on command line.\n"
         << endl;
    return 0;
  }

  string cmndfile = "";
  // Input command file.
  if(ip.hasOption("-c")) {
    cmndfile = ip.getOption("-c");
    if(cmndfile.find(".cmnd") == string::npos &&
        cmndfile.find(".dat") == string::npos) {
      cout << "Please provide a valid .cmnd file as "
      "argument to the -c option." << endl;
      return 1;
    }
  }
  else {
    cout << "You must provide a command file to produce output.\n"
            "Use option -h to show all command line options." << endl;
    return 1;
  }

  string cmndfile2 = "";
  // Optional secondary input command file.
  if(ip.hasOption("-c2")) {
    cmndfile2 = ip.getOption("-c2");
    if(cmndfile2.find(".cmnd") == string::npos &&
        cmndfile2.find(".dat") == string::npos) {
      cout << "Please provide a valid .cmnd file as argument "
      "to the -c2 option." << endl;
      return 1;
    }
  }

  string seed = "-1";
  // Optional seed from command line.
  if(ip.hasOption("-s")) seed = ip.getOption("-s");

  string out = "";
  // Set individual output prefix.
  if(ip.hasOption("-o")) out = ip.getOption("-o");

  int nev = -1;
  // Command line number of event, overrides the one set in input .cmnd file.
  if(ip.hasOption("-n")) nev = stoi(ip.getOption("-n"));

  // Catch the splash screen in a buffer.
  stringstream splashBuf;
  streambuf* sBuf = cout.rdbuf();
  cout.rdbuf(splashBuf.rdbuf());
  // The Pythia object.
  Pythia pythia;
  // Direct cout back.
  cout.rdbuf(sBuf);

  // Some extra parameters.
  pythia.settings.addFlag("Main:writeHepMC",false);
  pythia.settings.addFlag("Main:writeRoot",false);
  pythia.settings.addFlag("Main:runRivet",false);
  pythia.settings.addFlag("Main:rivetIgnoreBeams",false);
  pythia.settings.addFlag("Main:outputLog",false);
  pythia.settings.addWVec("Main:analyses",vector<string>());
  pythia.settings.addWVec("Main:preload",vector<string>());
  pythia.settings.addWord("Main:rivetRunName","");
  // Read input from external file.
  pythia.readFile(cmndfile);

  if(cmndfile2 != "") pythia.readFile(cmndfile2);

  // Set seed after reading input
  if(seed != "-1") {
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = "+seed);
  }

  // Read the extra parameters.
  int nEvent = pythia.mode("Main:numberOfEvents");;
  if(nev > -1) nEvent = nev;
  const bool hepmc = pythia.flag("Main:writeHepMC");
  const bool root = pythia.flag("Main:writeRoot");
  const bool runRivet = pythia.flag("Main:runRivet");
  const bool ignoreBeams = pythia.flag("Main:rivetIgnoreBeams");
  const bool doLog = pythia.flag("Main:outputLog");
  const string rivetrName = pythia.settings.word("Main:rivetRunName");
  const vector<string> rAnalyses = pythia.settings.wvec("Main:analyses");
  const vector<string> rPreload = pythia.settings.wvec("Main:preload");
  int nError = pythia.mode("Main:timesAllowErrors");
  bool countErrors = (nError > 0 ? true : false);
  // HepMC conversion object.
  HepMC::Pythia8ToHepMC ToHepMC;
  HepMC::IO_GenEvent* hepmcIO;
  if (hepmc)
    hepmcIO = new HepMC::IO_GenEvent((out == "" ? "pythia.hepmc"
      : out + ".hepmc"),ios::out);
  // Rivet initialization.
  Pythia8Rivet rivet(pythia,(out == "" ? "Rivet.yoda" : out + ".yoda"));
  rivet.ignoreBeams(ignoreBeams);
  for(int i = 0, N = rAnalyses.size(); i < N; ++i){
    string analysis = rAnalyses[i];
    size_t pos = analysis.find(":");
    // Simple case, no analysis parameters.
    if(pos == string::npos)
      rivet.addAnalysis(analysis);
    else {
      string an = analysis.substr(0,pos);
      analysis.erase(0, pos + 1);
      pos = analysis.find(":");
      string par = analysis.substr(0,pos);
      size_t pos2 = par.find("->");
      if (pos2 == string::npos){
         pythia.info.errorMsg("Error in main93: malformed"
          " parameter "+par);
      }
      string pKey = par.substr(0,pos2);
      string pVal = par.substr(pos2+2,par.length());
      rivet.addAnalysis(an+":"+pKey+"="+pVal);
    }
  }
  for(int i = 0, N = rPreload.size(); i < N; ++i)
    rivet.addPreload(rPreload[i]);
  rivet.addRunName(rivetrName);
  // Root initialization
  #ifdef USE_ROOT
  TFile* file;
  RootEvent* re;
  TTree* tree;
  #endif
  if (root) {
   // First test if root is available on system.
   #ifndef USE_ROOT
        cout << "Option Main::writeRoot = on requires a working,\n"
                "linked Root installation." << endl;
        return 1;
   #else
   string op = (out == "" ? "pythia.root" : out + ".root");
   file = TFile::Open(op.c_str(),"recreate" );
   re = new RootEvent();
   tree = new TTree("t","Pythia8 event tree");
   tree->Branch("events",&re);
   #endif
  }

  // Logfile initialization.
  ofstream logBuf;
  streambuf* oldCout;
  if(doLog) {
    oldCout = cout.rdbuf(logBuf.rdbuf());
    logBuf.open((out == "" ? "pythia.log" : out + ".log"));
  }
  // Option to trash the splash screen.
  ostream cnull(NULL);
  if(ip.hasOption("-l")) cnull << splashBuf.str();
  else cout << splashBuf.str();
  // Initialise Pythia.
  pythia.init();
  // Make a sanity check of initialized Rivet analyses
  if (!runRivet && rAnalyses.size() > 0 )
    pythia.info.errorMsg("Warning in main93: Rivet analyses initialized,"
                    "but runRivet set to off.");
  // Loop over events.
  for ( int iEvent = 0; iEvent < nEvent; ++iEvent ) {
    if ( !pythia.next() ) {
      if (countErrors && --nError < 0) {
        pythia.stat();
        cout << " \n *-------  PYTHIA STOPPED!  -----------------------*"
             << endl;
        cout << " | Event generation failed due to too many errors. |" << endl;
        cout << " *-------------------------------------------------*" << endl;
        return 1;
      }
      continue;
    }
    if (runRivet) rivet();
    if (hepmc) {
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      Info* pInfo = &pythia.info;
      if ( pInfo && pInfo->hiinfo ) {
      HepMC::HeavyIon ion;
      ion.set_Ncoll_hard(pInfo->hiinfo->nCollNDTot());
      ion.set_Ncoll(pInfo->hiinfo->nAbsProj() +
                    pInfo->hiinfo->nDiffProj() +
                    pInfo->hiinfo->nAbsTarg() +
                    pInfo->hiinfo->nDiffTarg() -
                    pInfo->hiinfo->nCollND() -
                    pInfo->hiinfo->nCollDD());
      ion.set_Npart_proj(pInfo->hiinfo->nAbsProj() +
                         pInfo->hiinfo->nDiffProj());
      ion.set_Npart_targ(pInfo->hiinfo->nAbsTarg() +
                         pInfo->hiinfo->nDiffTarg());
      ion.set_impact_parameter(pInfo->hiinfo->b());
      hepmcevt->set_heavy_ion(ion);
    }


      ToHepMC.fill_next_event( pythia, hepmcevt);
      (*hepmcIO) << hepmcevt;
      delete hepmcevt;
    }
    #ifdef USE_ROOT
    if (root) {
      // If we want to write a root file, the event must be skimmed here.
      vector<RootTrack> rts;
      for(int i = 0; i < pythia.event.size(); ++i) {
        RootTrack t;
        Particle& p = pythia.event[i];
        // Any particle cuts and further track definitions should
        // be implemented in the RootTrack class by the user.
        if (t.init(p)) rts.push_back(t);
      }
      bool fillTree = re->init(&pythia.info);
      re->tracks = rts;
      if(fillTree) tree->Fill();
    }
    #endif
    }
  if(hepmc) delete hepmcIO;

  pythia.stat();
  #ifdef USE_ROOT
  if (root) {
   tree->Print();
   tree->Write();
   delete file;
   delete re;
  }
  #endif
  // Put cout back in its place.
  if (doLog) cout.rdbuf(oldCout);
  return 0;
}
