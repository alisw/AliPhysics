// main38.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how Les Houches Event File version 3.0 information
// can be extracted from pythia.info.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator
  Pythia pythia;

  // Stick with default values, so do not bother with a separate file
  // for changes. However, do one change, to show readString in action.
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");

  // Initialize gzipped Les Houches Event File run.
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = wbj_lhef3.lhe");
  pythia.init();

  // Check if LHEF3 reader exists.
  bool hasRead = (pythia.info.LHEFversion() == 3);

  // Print header read from LHEF.
  cout << endl << "*****************************************************"
       << endl << "  PRINT HEADER OF LHE FILE" << endl;
  if (hasRead) cout << pythia.info.getHeaderBlock() << endl;
  cout << "*****************************************************" << endl;

  // Print extra LHEF v3 header information.
  cout << endl << "*****************************************************"
       << endl << "  RETRIEVE INITRWGT INFORMATION" << endl;
  if (hasRead) {
    LHAinitrwgt* initrwgt = pythia.info.initrwgt;
    cout << endl << "Attributes: " << endl;
    for ( map<string,string>::const_iterator
      it = initrwgt->attributes.begin();
      it != initrwgt->attributes.end(); ++it )
      cout << "  Attribute: " << it->first << ".  ." << it->second
           << endl;
    cout << "Contents: " << initrwgt->contents << endl;

    cout << endl << endl << "  RETRIEVE WEIGHTGROUP INFORMATION" << endl;
    map<string,LHAweightgroup> weightgroups = initrwgt->weightgroups;
    for ( map<string,LHAweightgroup>::const_iterator
      it = weightgroups.begin();
      it != weightgroups.end(); ++it ) {
      cout << "  Weightgroup: " << it->first << endl;
      for ( map<string,string>::const_iterator
        it2 = it->second.attributes.begin();
        it2 != it->second.attributes.end(); ++it2 ) {
        cout << "  Attribute: " << it2->first << ".  ." << it2->second
             << endl;
      }
      cout << endl << "  Retrieve weight information:" << endl;
      map<string,LHAweight> weights = it->second.weights;
      for ( map<string,LHAweight>::const_iterator
        it2 = weights.begin(); it2 != weights.end(); ++it2 ) {
        cout << "  Weight: " << it2->first << endl;
        for ( map<string,string>::const_iterator
          it3 = it2->second.attributes.begin();
          it3 != it2->second.attributes.end(); ++it3 ) {
          cout << "  Attribute: " << it3->first << ".  ." << it3->second
               << endl;
        }
        cout << "  Contents: " << it2->second.contents << endl;
      }
    }
  }
  cout << "*****************************************************" << endl;

  // Print extra LHEF v3 initialization information.
  cout << endl << "*****************************************************"
       << endl << "  RETRIEVE GENERATOR INFORMATION" << endl;
  unsigned int ngen = 0;
  if (hasRead) ngen = pythia.info.generators->size();
  cout << "Number of generator tags " << ngen << endl;
  for (unsigned int igen = 0; igen < ngen; ++igen) {
    cout << ".  generator tag ." << igen << ". : ."
         << (*pythia.info.generators)[igen].contents << ".\n"
         << ".  name          ."
         << (*pythia.info.generators)[igen].attributes["name"]
         << ".  version       ."
         << (*pythia.info.generators)[igen].attributes["version"]
         << ".  random        ."
         << (*pythia.info.generators)[igen].attributes["random"]
         << "." << endl;
  }
  cout << "*****************************************************" << endl;

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {

    cout << endl << "Read event # " << iEvent << endl;

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;
      else continue;
    }

    // Print process-level event.
    pythia.process.list();

    // Print extra LHEF v3 event information.
    cout << endl << "Print event # " << iEvent << endl;

    cout << endl << "*****************************************************"
         << endl << "  RETRIEVE EVENT INFORMATION" << endl;
    if (hasRead) {
      map <string,string> eventAttributes =
        *(pythia.info.eventAttributes);
      cout << endl << "Attributes:" << endl;
      for ( map<string,string>::const_iterator
        it = eventAttributes.begin(); it != eventAttributes.end(); ++it )
          cout << "  Attribute\t" << it->first << " = " << it->second
               << endl;
    }
    cout << "*****************************************************" << endl;

    // Print extra LHEF v3 event weight information.
    cout << endl << "*****************************************************"
         << endl << "  RETRIEVE WEIGHTS (DETAILED FORMAT) INFORMATION" << endl;
    if (hasRead) {
      int nwgt = pythia.info.rwgt->size();
      cout << "Number of wgt tags " << nwgt << endl;
      map<string,double> weights = *(pythia.info.weights_detailed);
      for ( map<string,double>::const_iterator it = weights.begin();
        it != weights.end(); ++it ) {
        cout << it->first << " " << it->second << " " << endl;
      }
      cout << "Second option, from ordered vector:" << endl;
      for ( int i = 0;
        i < int(pythia.info.weights_detailed_vector.size());
        ++i ) {
        cout << "Weight " << i << " : "
             << pythia.info.weights_detailed_vector[i] << endl;
      }
    }
    cout << "*****************************************************" << endl;

    // Print extra LHEF v3 compressed format event weight information.
    cout << endl << "*****************************************************"
         << endl << "  RETRIEVE WEIGHTS (COMPRESSED FORMAT) INFORMATION"
         << endl;
    if (hasRead) {
      int nweights = pythia.info.weights->size();
      cout << "Number of weights (only one tag!) " << nweights << endl;
      LHAweights* weights = pythia.info.weights;
      cout << endl << "Attributes:" << endl;
      for ( map<string,string>::const_iterator
        it = weights->attributes.begin();
        it != weights->attributes.end(); ++it )
          cout << "  Attribute\t" << it->first << " = " << it->second << endl;
      cout << "Weights:" << endl;
      for ( int i = 0; i < int(weights->weights.size()); ++i)
        cout << "  Weight\t" << weights->weights[i] << endl;
      cout << endl << "Contents:" << weights->contents << endl;
      cout << "Second option, from stored vector:" << endl;
      for ( int i = 0;
        i < int(pythia.info.weights_compressed->size());
        ++i ) {
        cout << "Weight\t" << i << " : "
             << (*pythia.info.weights_compressed)[i] << endl;
      }
    }
    cout << "*****************************************************" << endl;

    // Print extra LHEF v3 scale information.
    cout << endl << "*****************************************************"
         << endl << "  RETRIEVE SCALES INFORMATION" << endl;
    if (hasRead) {
      LHAscales* scales = pythia.info.scales;
      cout << endl << "Attributes:" << endl;
      for ( map<string,double>::const_iterator
        it = scales->attributes.begin(); it != scales->attributes.end(); ++it )
          cout << "  Attribute\t" << it->first << " = " << it->second << endl;
      cout << "Contents:" << scales->contents << endl;
    }
    cout << "*****************************************************" << endl;

    // Print event comments.
    cout << endl << "*****************************************************"
         << endl << "  RETRIEVE EVENT COMMENTS" << endl;
    if (hasRead) cout << pythia.info.getEventComments() << endl;
    cout << "*****************************************************" << endl;


  // End of event loop.
  }

  // Give statistics.
  pythia.stat();

  // Done.
  return 0;
}
