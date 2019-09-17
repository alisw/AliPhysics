// main33.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Philip Ilten.

// An example where the HVQ POWHEGBOX matrix element binary is
// interfaced directly with PYTHIA. For this example to run correctly
// PYTHIA must be configured with the
// --with-powheg-bin=<path to directory containing only POWHEG binaries>
// option. This will build plugin libraries of the name
// libpythia8powheg<binary name>.so in the library directory.
// For these plugin libraries to build correctly, special compiler flags
// must have been used when building the POWHEGBOX binaries. These are
// "-rdynamic -fPIE -fPIC -pie". The following SED command will correctly
// insert them into the relevant POWHEGBOX Makefile:
//     sed -i "s/F77= gfortran/F77= gfortran -rdynamic -fPIE -fPIC -pie/g"
//     Makefile
// For this specific example the library libpythia8powheghvq.so must
// have been built.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegProcs.h"

using namespace Pythia8;

int main() {

  Pythia pythia;

  // The constructor PowhegProcs(process name, pythia, run directory,
  // PDF filename, use Pythia random) where run directory by default
  // is ./powhegrun, the PDF filename is empty, and using Pythia
  // random is true. If using native PDFs rather than LHAPDF, the PDF
  // filename is the full name of the PDF file to copy to the run
  // directory.
  PowhegProcs procs(&pythia, "hvq");

  // The PowhegProcs class automatically sets the Pythia user hooks to
  // an instance of PowhegHooks. However, this can be modified to a
  // user chosen set of hooks (or null).
  // pythia->setUserHooksPtr(userHooksPtr);

  // Pythia and the POWHEG user hooks must still be configured, here
  // this is done via main33.cmnd. These settings are sensible
  // defaults, but Powheg:nFinal is dependent upon the POWHEG matrix
  // element being used and so must be changed as appropriate.
  pythia.readFile("main33.cmnd");

  // The commands readFile and readString are used to configure the
  // POWHEG matrix element. If a setting is repeated a warning is
  // issued and the most recent setting is used.
  procs.readFile("main33.pwhg");

  // This init call must be made before PYTHIA is initialized. It
  // copies the POWHEG input and PDF file to the POWHEG run directory.
  procs.init();

  // Initialize Pythia, based on the specified settings.
  pythia.init();

  // Run PYTHIA. The random numbers are taken from the associated
  // PYTHIA random number generator.
  for (int iEvent = 0; iEvent < 100; ++iEvent) pythia.next();

  // End of run.
  pythia.stat();
  return 0;
}
