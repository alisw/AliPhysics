// PowhegProcs.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Philip Ilten, May 2015.

#ifndef Pythia8_PowhegProcs_H
#define Pythia8_PowhegProcs_H

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"

namespace Pythia8 {

//==========================================================================

// A class to generate events with hard processes from POWHEGBOX
// matrix elements. See http://powhegbox.mib.infn.it/ for further
// details on POWHEGBOX.

// WARNING: If one wishes to use LHAPDF with both POWHEGBOX and
// Pythia, only LHAPDF 6 configured with the option
// "--with-lhapdf6-plugin=LHAPDF6.h" for Pythia should be used. If
// not, and differing PDF sets are used between POWHEGBOX and Pythia,
// POWHEGBOX will not re-initialize the PDF set and consequently will
// use the PDF set last used by Pythia.

class PowhegProcs {

public:

  // Constructor.
  PowhegProcs(Pythia *pythiaPtrIn, string procIn, string dirIn = "powhegrun",
    string pdfIn = "", bool random = true);

  // Destructor.
  ~PowhegProcs();

  // Read a POWHEG settings string.
  bool readString(string line);

  // Read a POWHEG settings file.
  bool readFile(string name);

  // Write out the input for POWHEG.
  bool init();

  // The POWHEG LHAup pointer.
  LHAup *lhaup;

private:

  // Typedefs of the hooks used to access the plugin.
  typedef LHAup* NewLHAupPowheg(Pythia*);
  typedef void DeleteLHAupPowheg(LHAup*);

  // The POWHEG process name, run directory, and PDF file (if not LHAPDF).
  string proc, dir, pdf;

  // The map of POWHEG settings.
  map<string, string> settings;

  // The associated PYTHIA pointer.
  Pythia *pythia;

  // The POWHEG plugin library.
  void *lib;

  // The POWHEG hooks.
  PowhegHooks hooks;

};

//--------------------------------------------------------------------------

// Constructor.

// pythiaPtrIn: The PYTHIA object the plugin will use for settings and
// random numbers.

// procIn: the process name. An attempt is made to load the plugin
// library libpythia8powheg<procIn>.so.

// dirIn: The directory where the POWHEG matrix element will be
// run. This is needed if two instances are to be run concurrently
// since the matrix element generates a large number of files.

// pdfIn: The full path and name of the PDF file to use, if not using
// LHAPDF. This file is copied to the run directory via the init()
// method.

// random: Flag to use the Pythia random number generator with
// POWHEGBOX. If true, the POWHEGBOX random number block is
// initialized for each event from the Pythia random number
// generator. If false, the default POWHEGBOX random number generation
// is performed. Note that the initialization is always performed
// using the POWHEGBOX random number generation (which can be modified
// via the POWHEGBOX configuration).

PowhegProcs::PowhegProcs(Pythia *pythiaPtrIn, string procIn, string dirIn,
  string pdfIn, bool random) : lhaup(0), proc(procIn), dir(dirIn), pdf(pdfIn),
  pythia(pythiaPtrIn), lib(0) {

  if (!pythia) return;
  NewLHAupPowheg *sym(0);
  const char* error(0);

  // Load the library.
  lib = dlopen(("libpythia8powheg" + proc + ".so").c_str(), RTLD_LAZY);
  error = dlerror();
  if (error) {
    pythia->info.errorMsg("Error from PowhegProcs::PowhegProcs: "
      + string(error));
    return;
  }
  dlerror();

  // Load the LHAup pointer.
  sym = (NewLHAupPowheg*)dlsym(lib, "newLHAupPowheg");
  error = dlerror();
  if (error) {
    pythia->info.errorMsg("Error from PowhegProcs::PowhegProcs: "
      + string(error));
  }
  dlerror();
  pythia->settings.addWord("POWHEG:dir", dir);
  pythia->settings.addFlag("POWHEG:pythiaRandom", random);
  if (sym) lhaup = sym(pythia);

  // Configure PYTHIA.
  pythia->setLHAupPtr(lhaup);
  pythia->setUserHooksPtr(&hooks);

}

//--------------------------------------------------------------------------

// Destructor.

PowhegProcs::~PowhegProcs() {

  // Delete the LHAup pointer.
  if (lhaup && lib) {
    DeleteLHAupPowheg *sym(0);
    sym = (DeleteLHAupPowheg*)dlsym(lib, "deleteLHAupPowheg");
    if (sym) sym(lhaup);
  }

  // Unload the library.
  if (lib) {dlclose(lib); dlerror();}

}

//--------------------------------------------------------------------------

// Read a POWHEG settings string. If a setting is repeated a warning
// is printed but the most recent setting is used.

bool PowhegProcs::readString(string line) {

  // Copy string without initial and trailing blanks.
  if (!pythia) return false;
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return true;
  int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");
  int lastChar  = line.find_last_not_of(" \n\t\v\b\r\f\a");
  line = line.substr(firstChar, lastChar + 1 - firstChar);

  // Find the key.
  firstChar = line.find_first_of("  \t\f\v\n\r");
  string key = toLower( line.substr(0, firstChar), false);

  // Add the setting.
  if (key.size() > 0
    && key.find_first_of("abcdedfghijklmnopqrtsuvwxyz") == 0) {
    map<string, string>::iterator setting = settings.find(key);
    if (setting != settings.end()) {
      pythia->info.errorMsg("Warning from PowhegProcs::readString: replacing "
        "previous POWHEG setting for " + key + ".");
      setting->second = line;
    } else settings[key] = line;
  }
  return true;

}

//--------------------------------------------------------------------------

// Read a POWHEG settings file.

bool PowhegProcs::readFile(string name) {

  fstream config(name.c_str(), ios::in); string line;
  while (getline(config, line, '\n')) readString(line);
  config.close();
  return true;

}

//--------------------------------------------------------------------------

// Write the input for POWHEG.

bool PowhegProcs::init() {

  // Copy over the PDF file if needed.
  if (pdf != "") {
    fstream pdfin(pdf.c_str(), ios::in | ios::binary);
    fstream pdfout((dir + "/" + pdf.substr(0, pdf.find_last_of("/"))).c_str(),
      ios::out | ios::binary);
    pdfout << pdfin.rdbuf();
    pdfin.close();
    pdfout.close();
  }

  // Copy the settings to the configuration file.
  fstream config((dir + "/" + "powheg.input").c_str(), ios::out);
  for (map<string, string>::iterator setting = settings.begin();
    setting != settings.end(); ++setting) config << setting->second << "\n";
  config.close();
  return true;

}

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_PowhegProcs_H
