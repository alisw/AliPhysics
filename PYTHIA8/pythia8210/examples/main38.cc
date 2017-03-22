// main38.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
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
  // pythia.readString("PartonLevel:ISR = off");
  // pythia.readString("PartonLevel:FSR = off");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");
  //pythia.readString("PDF:pSet = LHAPDF5:cteq6m.LHpdf");

  // Initialize gzipped Les Houches Event File run.
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = wbj_lhef3.lhe");
  pythia.init();

  // Print extra LHEF v3 initialization information.
  cout << endl << "*****************************************************"
       << endl << "  RETRIEVE GENERATOR INFORMATION" << endl;
  unsigned int ngen = pythia.info.getGeneratorSize();
  cout << "Number of generator tags " << ngen << endl;
  for (unsigned int igen = 0; igen < ngen; ++igen) {
    cout << ".  generator tag ." << igen << ". : ."
         << pythia.info.getGeneratorValue(igen) << ".\n"
         << ".  name          ."
         << pythia.info.getGeneratorAttribute(igen,"name",igen)
         << ".  version       ."
         << pythia.info.getGeneratorAttribute(igen,"version",igen)
         << ".  random        ."
         << pythia.info.getGeneratorAttribute(igen,"random",igen)
         << "." << endl;
    cout << ".  generator tag ." << igen << ". : ."
         << pythia.info.getGeneratorValue(igen) << ".\n"
         << ".  name          ."
         << pythia.info.getGeneratorAttribute(igen,"name",true)
         << ".  version       ."
         << pythia.info.getGeneratorAttribute(igen,"version",true)
         << ".  random        ."
         << pythia.info.getGeneratorAttribute(igen,"random",true)
         << "." << endl;
  }

  // Get number of event weights.
  int ninitrwgt = pythia.info.getInitrwgtSize();

  // Initialise as many histograms as there are event weights.
  vector<Hist> pTw;
  for (int iHist = 0; iHist < ninitrwgt; ++iHist)
    pTw.push_back( Hist("pT W",50,0.,200.) );

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
    if (pythia.info.eventAttributes) {
      for ( map<string,string>::const_iterator
        it = pythia.info.eventAttributes->begin();
        it != pythia.info.eventAttributes->end(); ++it )
          cout << ".  attributes    ." << it->first << ".  ." << it->second
               << ".  .";
      cout << "." << endl;
      cout << pythia.info.getEventAttribute("npLO") << ".  ."
           << pythia.info.getEventAttribute("npNLO") << endl;
      cout << pythia.info.getEventAttribute("npLO",true) << ".  ."
           << pythia.info.getEventAttribute("npNLO",true) << endl;
    }

    // Print extra LHEF v3 event weight information.
    cout << endl << "*****************************************************"
         << endl << "  RETRIEVE WEIGHTS (DETAILED FORMAT) INFORMATION" << endl;
    unsigned int nwgt = pythia.info.getWeightsDetailedSize();
    cout << "Number of wgt tags " << nwgt << endl;
    for (unsigned int iwgt = 0; iwgt < nwgt; ++iwgt) {
      string key;
      ostringstream convert;
      convert << iwgt + 1001;
      key = convert.str();
      cout << ".  wgt tag       ." << iwgt << ". : ."
           << pythia.info.getWeightsDetailedValue(key) << ".\n"
           << ".  id            ."
           << pythia.info.getWeightsDetailedAttribute(key,"id")
           << ".  version       ."
           << pythia.info.getWeightsDetailedAttribute(key,"version")
           << "." << endl;
    }

    // Print extra LHEF v3 event weight information directly from iterator.
    cout << "wgt tags directly from iterator" << endl;
    unsigned int nwgts = 0;
    if (pythia.info.rwgt) nwgts = pythia.info.rwgt->wgts.size();
    cout << "Number of weight tags " << nwgts << endl;
    for (unsigned int iwgts = 0; iwgts < nwgts; ++iwgts) {
      string key;
      ostringstream convert;
      convert << iwgts+1001;
      key = convert.str();
      cout << ".  wgt tag       ." << iwgts << ". : ."
           << pythia.info.rwgt->wgts[key].contents << ".\n"
           << ".  id            ." << pythia.info.rwgt->wgts[key].id;
      for ( std::map<std::string,std::string>::const_iterator
        it = pythia.info.rwgt->wgts[key].attributes.begin();
        it != pythia.info.rwgt->wgts[key].attributes.end(); ++it )
          cout << ".  attributes    ." << it->first << ".  ." << it->second
               << ".  .";
      cout << "." << endl;
    }

    // Print extra LHEF v3 compressed format event weight information.
    cout << endl << "*****************************************************"
         << endl << "  RETRIEVE WEIGHTS (COMPRESSED FORMAT) INFORMATION"
         << endl;
    unsigned int nweights = pythia.info.getWeightsCompressedSize();
    cout << "Number of weights (only one tag!) " << nweights << endl;
    for (unsigned int iweights = 0; iweights < nweights; ++iweights)
      cout << ".  weight        ." << iweights << ". : ."
           << pythia.info.getWeightsCompressedValue(iweights) << ".\n"
           << "." << endl;

    // Print extra LHEF v3 weight tags directly from iterator.
    cout << "weight tags directly from iterator" << endl;
    if (pythia.info.weights)
    for ( std::map<std::string,std::string>::const_iterator
      it = pythia.info.weights->attributes.begin();
      it != pythia.info.weights->attributes.end(); ++it )
        cout << ".  attributes    ." << it->first << ".  ." << it->second
             << ".  .";
    cout << "." << endl;

    // Print extra LHEF v3 scale information.
    cout << endl << "*****************************************************"
         << endl << "  RETRIEVE SCALES INFORMATION" << endl;
    cout << ".  ." << pythia.info.getScalesValue() << endl;
    cout << ".  ." << pythia.info.getScalesAttribute("muf") << endl;
    cout << ".  ." << pythia.info.getScalesAttribute("mur") << endl;
    cout << ".  ." << pythia.info.getScalesAttribute("mups") << endl;
    cout << ".  ." << pythia.info.getScalesAttribute("SCALUP") << endl;
    if(pythia.info.scales)
    for ( std::map<std::string,double>::const_iterator
      it = pythia.info.scales->attributes.begin();
      it != pythia.info.scales->attributes.end(); ++it )
        cout << ".  attributes    ." << it->first << ".  ." << it->second
             << ".  .";
    cout << "." << endl;

    // Find the final copy of the W.
    int iW = 0;
    for (int i = pythia.event.size()-1; i > 0; --i)
      if (pythia.event[i].idAbs() == 24) { iW = i; break;}
    double pT = pythia.event[iW].pT();

    // Loop over the event weights in the detailed format and histogram.
    nwgt = pythia.info.getWeightsDetailedSize();
    for (unsigned int iwgt = 0; iwgt < nwgt; ++iwgt) {
      string key;
      ostringstream convert;
      convert << iwgt + 1001;
      key = convert.str();
      double w = pythia.info.getWeightsDetailedValue(key);
      pTw[iwgt].fill( max(pT,0.5), w );
    }

  // End of event loop.
  }

  // Give statistics.
  pythia.stat();

  // Print histograms.
  ofstream write;
  stringstream suffix;
  for (int iHist = 0; iHist < ninitrwgt; ++iHist) {
    suffix << iHist << ".dat";
    // Write histograms to file.
    write.open( (char*)("PTW_" + suffix.str()).c_str());
    pTw[iHist].table(write);
    suffix.str("");
    write.close();
  }

  // Done.
  return 0;
}
