// MergingHooks.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel.
// Function definitions (not found in the header) for the HardProcess and
// MergingHooks classes.

#include "Pythia8/MergingHooks.h"

namespace Pythia8 {

//==========================================================================

// The HardProcess class.

//--------------------------------------------------------------------------

// Declaration of hard process class
// This class holds information on the desired hard 2->2 process to be merged
// This class is a container class for History class use.

// Initialisation on the process string

void HardProcess::initOnProcess( string process, ParticleData* particleData) {
  state.init("(hard process)", particleData);
  translateProcessString(process);
}

//--------------------------------------------------------------------------

// Initialisation on the path to LHE file

void HardProcess::initOnLHEF( string LHEfile, ParticleData* particleData) {
  state.init("(hard process)", particleData);
  translateLHEFString(LHEfile);
}

//--------------------------------------------------------------------------

// Function to access the LHE file and read relevant information.
// The merging scale will be read from the +1 jet sample, called
//   LHEpath_1.lhe
// while the hard process will be read from
//   LHEpath_0.lhe
// Currently, only read from MadEvent- or Sherpa-generated LHE files
// is automatic, else the user is asked to supply the necessary
// information.

void HardProcess::translateLHEFString( string LHEpath){

  // Open path to LHEF and extract merging scale
  ifstream infile;
  infile.open( (char*)( LHEpath +"_0.lhe").c_str());

  // Check with ME generator has been used
  int iLine =0;
  int nLinesMax = 200;
  string lineGenerator;
  while( iLine < nLinesMax
        && lineGenerator.find("SHERPA", 0) == string::npos
        && lineGenerator.find("POWHEG-BOX", 0) == string::npos
        && lineGenerator.find("Pythia8", 0) == string::npos
        && lineGenerator.find("MadGraph", 0) == string::npos){
    iLine++;
    lineGenerator = " ";
    getline(infile,lineGenerator);
  }
  infile.close();

  vector <int> incom;
  vector <int> inter;
  vector <int> outgo;
  // Particle identifiers, ordered in such a way that e.g. the "u"
  // in a mu is not mistaken for an u quark
  int inParticleNumbers[] = {
        // Leptons
        -11,11,-12,12,-13,13,-14,14,-15,15,-16,16,
        // Jet container
        2212,2212,0,0,0,0,
        // Quarks
        -1,1,-2,2,-3,3,-4,4,-5,5,-6,6};

  string inParticleNamesSH[] = {
        // Leptons
        "-11","11","-12","12","-13","13","-14","14","-15","15","-16","16",
        // Jet container
        "-93","93","-90","90","-91","91",
        // Quarks
        "-1","1","-2","2","-3","3","-4","4","-5","5","-6","6"};
  string inParticleNamesMG[] =  {
        // Leptons
        "e+","e-","ve~","ve","mu+","mu-","vm~","vm","ta+","ta-","vt~","vt",
        // Jet container
        "p~","p","l+","l-","vl~","vl",
        // Quarks
        "d~","d","u~","u","s~","s","c~","c","b~","b","t~","t"};

  // Declare intermediate particle identifiers
  int interParticleNumbers[] = {
         // Electroweak gauge bosons
         22,23,-24,24,25,2400,
         // Top quarks
         -6,6,
         // Dummy index as back-up
         0,
         // All squarks
        -1000001,1000001,-1000002,1000002,-1000003,1000003,-1000004,1000004,
        -1000005,1000005,-1000006,1000006,-2000001,2000001,-2000002,2000002,
        -2000003,2000003,-2000004,2000004,-2000005,2000005,-2000006,2000006};
  // Declare names of intermediate particles
  string interParticleNamesMG[] = {
        // Electroweak gauge bosons
        "a","z","w-","w+","h","W",
         // Top quarks
         "t~","t",
        // Dummy index as back-up
        "xx",
         // All squarks
         "dl~","dl","ul~","ul","sl~","sl","cl~","cl","b1~","b1","t1~","t1",
         "dr~","dr","ur~","ur","sr~","sr","cr~","cr","b2~","b2","t2~","t2"};

  // Declare final state particle identifiers
  int outParticleNumbers[] = {
        // Leptons
        -11,11,-12,12,-13,13,-14,14,-15,15,-16,16,
        // Jet container and lepton containers
        2212,2212,0,0,0,0,1200,1100,5000,
        // Quarks
        -1,1,-2,2,-3,3,-4,4,-5,5,-6,6,
        // SM uncoloured bosons
        22,23,-24,24,25,2400,
        // Neutralino in SUSY
        1000022,
        // All squarks
        -1000001,1000001,-1000002,1000002,-1000003,1000003,-1000004,1000004,
        -1000005,1000005,-1000006,1000006,-2000001,2000001,-2000002,2000002,
        -2000003,2000003,-2000004,2000004,-2000005,2000005,-2000006,2000006};
  // Declare names of final state particles
  string outParticleNamesMG[] =  {
        // Leptons
        "e+","e-","ve~","ve","mu+","mu-","vm~","vm","ta+","ta-","vt~","vt",
        // Jet container and lepton containers
        "j~","j","l+","l-","vl~","vl","NEUTRINOS","LEPTONS","BQUARKS",
        // Quarks
        "d~","d","u~","u","s~","s","c~","c","b~","b","t~","t",
        // SM uncoloured bosons
        "a","z","w-","w+","h","W",
        // Neutralino in SUSY
        "n1",
        // All squarks
        "dl~","dl","ul~","ul","sl~","sl","cl~","cl","b1~","b1","t1~","t1",
        "dr~","dr","ur~","ur","sr~","sr","cr~","cr","b2~","b2","t2~","t2"};

  string outParticleNamesSH[] = {
        // Leptons
        "-11","11","-12","12","-13","13","-14","14","-15","15","-16","16",
        // Jet container and lepton containers
        "-93","93","-90","90","-91","91","0","0","0",
        // Quarks
        "-1","1","-2","2","-3","3","-4","4","-5","5","-6","6",
        // SM uncoloured bosons
        "22","23","-24","24","25","0",
        // Neutralino in SUSY
        "1000022",
        // All squarks
        "-1000001","1000001","-1000002","1000002","-1000003","1000003",
                   "-1000004","1000004",
        "-1000005","1000005","-1000006","1000006","-2000001","2000001",
                   "-2000002","2000002",
        "-2000003","2000003","-2000004","2000004","-2000005","2000005",
                   "-2000006","2000006"};

  // Declare size of particle name arrays
  int nIn   = 30;
  int nInt  = 33;
  int nOut  = 64;

  // Save type of the generator, in order to be able to extract
  // the tms definition
  int meGenType = (lineGenerator.find("MadGraph", 0) != string::npos) ? -1
              : (lineGenerator.find("SHERPA", 0) != string::npos)     ? -2
              : (lineGenerator.find("POWHEG-BOX", 0) != string::npos) ? -3
              : (lineGenerator.find("Pythia8", 0) != string::npos)    ? -4
              : 0;

  if (meGenType == -2){
    // Now read merging scale
    // Open path to LHEF and extract merging scale
    infile.open( (char*)( LHEpath +"_1.lhe").c_str());
    string lineTMS;
    while(lineTMS.find("NJetFinder ", 0) == string::npos){
      lineTMS = " ";
      getline(infile,lineTMS);
    }
    infile.close();
    lineTMS = lineTMS.substr(0,lineTMS.find(" 0.0 ",0));
    lineTMS = lineTMS.substr(lineTMS.find(" ", 0)+3,lineTMS.size());
    // Remove whitespaces
    while(lineTMS.find(" ", 0) != string::npos)
      lineTMS.erase(lineTMS.begin()+lineTMS.find(" ",0));
    // Replace d with e
    if ( lineTMS.find("d", 0) != string::npos)
      lineTMS.replace(lineTMS.find("d", 0),1,1,'e');
    tms = atof((char*)lineTMS.c_str());

    // Now read hard process
    // Open path to LHEF and extract hard process
    infile.open( (char*)( LHEpath +"_0.lhe").c_str());
    string line;
    while(line.find("Process", 0) == string::npos){
      line = " ";
      getline(infile,line);
    }
    infile.close();
    line = line.substr(line.find(" ",0),line.size());

    // Cut string into incoming and outgoing pieces
    vector <string> pieces;
    pieces.push_back( line.substr(0,line.find("->", 0)) );
    // Do not count additional final jets
    int end = (line.find("{", 0) != string::npos) ? line.find("{", 0)-2
            : line.size();
    pieces.push_back( line.substr(line.find("->", 0)+2, end) );

    // Get incoming particles
    for(int i=0; i < nIn; ++i) {
      for(int n = pieces[0].find(inParticleNamesSH[i], 0);
             n != int(string::npos);
             n = pieces[0].find(inParticleNamesSH[i], n)) {
        incom.push_back(inParticleNumbers[i]);
        pieces[0].erase(pieces[0].begin()+n,
                        pieces[0].begin()+n+inParticleNamesSH[i].size());
        n=0;
      }
    }
    // Get intermediate particles
    // If intermediates are still empty, fill intermediate with default value
    inter.push_back(0);
    // Get final particles
    for(int i=0; i < nOut; ++i) {
      for(int n = pieces[1].find(outParticleNamesSH[i], 0);
             n != int(string::npos);
             n = pieces[1].find(outParticleNamesSH[i], n)) {
        outgo.push_back(outParticleNumbers[i]);
        pieces[1].erase(pieces[1].begin()+n,
                        pieces[1].begin()+n+outParticleNamesSH[i].size());
        n=0;
      }
    }

  } else if (meGenType == -1 || meGenType == -3 || meGenType == -4){

    // Now read merging scale
    string lineTMS;

    if (meGenType == -1) {
      // Open path to LHEF and extract merging scale
      infile.open( (char*)( LHEpath +"_1.lhe").c_str());
      while(lineTMS.find("ktdurham", 0) == string::npos){
        lineTMS = " ";
        getline(infile,lineTMS);
      }
      lineTMS = lineTMS.substr(0,lineTMS.find("=",0));
      infile.close();
    } else {
      lineTMS = "30.";
    }

    // Remove whitespaces
    while(lineTMS.find(" ", 0) != string::npos)
      lineTMS.erase(lineTMS.begin()+lineTMS.find(" ",0));
    // Replace d with e
    if ( lineTMS.find("d", 0) != string::npos)
      lineTMS.replace(lineTMS.find("d", 0),1,1,'e');
    tms = atof((char*)lineTMS.c_str());

    // Now read hard process
    // Open path to LHEF and extract hard process
    infile.open( (char*)( LHEpath +"_0.lhe").c_str());
    string line;
    while(line.find("@1", 0) == string::npos){
      line = " ";
      getline(infile,line);
    }
    infile.close();
    line = line.substr(0,line.find("@",0));

    // Count number of resonances
    int appearances = 0;
    for(int n = line.find("(", 0); n != int(string::npos);
            n = line.find("(", n)) {
      appearances++;
      n++;
    }

    // Cut string in incoming, resonance+decay and outgoing pieces
    vector <string> pieces;
    for(int i =0; i < appearances;++i) {
      int n = line.find("(", 0);
      pieces.push_back(line.substr(0,n));
      line = line.substr(n+1,line.size());
    }
    // Cut last resonance from rest
    if (appearances > 0) {
      pieces.push_back( line.substr(0,line.find(")",0)) );
      pieces.push_back( line.substr(line.find(")",0)+1,line.size()) );
    }

    // If the string was not cut into pieces, i.e. no resonance was
    // required, cut string using '>' as delimiter
    if (pieces.empty() ){
      appearances = 0;
      for(int n = line.find(">", 0); n != int(string::npos);
              n = line.find(">", n)) {
        appearances++;
        n++;
      }

      // Cut string in incoming and outgoing pieces
      for(int i =0; i < appearances;++i) {
        int n = line.find(">", 0);
        pieces.push_back(line.substr(0,n));
        line = line.substr(n+1,line.size());
      }

      if (appearances == 1) pieces.push_back(line);
      if (appearances > 1) {
        pieces.push_back( line.substr(0,line.find(">",0)) );
        pieces.push_back( line.substr(line.find(">",0)+1,line.size()) );
      }
    }

    // Get incoming particles
    for(int i=0; i < nIn; ++i) {
      for(int n = pieces[0].find(inParticleNamesMG[i], 0);
             n != int(string::npos);
             n = pieces[0].find(inParticleNamesMG[i], n)) {
        incom.push_back(inParticleNumbers[i]);
        pieces[0].erase(pieces[0].begin()+n,
                        pieces[0].begin()+n+inParticleNamesMG[i].size());
        n=0;
      }
    }

    // Check intermediate resonances and decay products
    for(int i =1; i < int(pieces.size()); ++i){
      // Seperate strings into intermediate and outgoing, if not already done
      int k = pieces[i].find(">", 0);

      string intermediate = (pieces[i].find(">", 0) != string::npos) ?
                             pieces[i].substr(0,k) : "";
      string outgoing = (pieces[i].find(">", 0) != string::npos) ?
                         pieces[i].substr(k+1,pieces[i].size()) : pieces[i];

      // Get intermediate particles
      for(int j=0; j < nInt; ++j) {
        for(int n = intermediate.find(interParticleNamesMG[j], 0);
               n != int(string::npos);
               n = intermediate.find(interParticleNamesMG[j], n)) {
          inter.push_back(interParticleNumbers[j]);
          intermediate.erase(intermediate.begin()+n,
                      intermediate.begin()+n+interParticleNamesMG[j].size());
          n=0;
        }
      }

      // Get outgoing particles
      for(int j=0; j < nOut; ++j) {
        for(int n = outgoing.find(outParticleNamesMG[j], 0);
               n != int(string::npos);
               n = outgoing.find(outParticleNamesMG[j], n)) {
          outgo.push_back(outParticleNumbers[j]);
          outgoing.erase(outgoing.begin()+n,
                         outgoing.begin()+n+outParticleNamesMG[j].size());
          n=0;
        }
      }

      // For arbitrary or non-existing intermediate, remember zero for each
      // two outgoing particles, without bosons.
      if (inter.empty()) {

        // For final state bosons, bookkeep the final state boson as
        // intermediate as well.
        int nBosons = 0;
        for(int l=0; l < int(outgo.size()); ++l)
          if ( (abs(outgo[l]) > 20 && abs(outgo[l]) <= 25) || outgo[l] == 2400)
            nBosons++;

        int nZeros = (outgo.size() - nBosons)/2;
        for(int l=0; l < nZeros; ++l)
          inter.push_back(0);
      }

      // For final state bosons, bookkeep the final state boson as
      // intermediate as well.
      for(int l=0; l < int(outgo.size()); ++l)
        if ( (abs(outgo[l]) > 20 && abs(outgo[l]) <= 25) || outgo[l] == 2400)
          inter.push_back(outgo[l]);

    }

  } else {

    cout << "Reading of tms and hard process information from LHEF currently"
         << " only automated for MadEvent- or SHERPA-produced LHEF" << endl;
    int tempInt      = 0;
    cout << "Use default process pp -> e+ve + jets? (0:no / 1:yes): ";
    cin >> tempInt;
    cout << endl;

    if (tempInt == 0){
      tempInt = 0;
      double tempDouble  = 0.0;
      cout << "Please specify merging scale (kT Durham, in GeV): ";
      cin >> tempDouble;
      tms = tempDouble;
      cout << endl;
      cout << "Please specify first incoming particle ";
      cout << "(p+/p- = 2212, e- = 11, e+ = -11): ";
      cin >> tempInt;
      incom.push_back(tempInt);
      tempInt = 0;
      cout << endl;
      cout << "Please specify second incoming particle ";
      cout << "(p+/p- = 2212, e- = 11, e+ = -11): ";
      cin >> tempInt;
      incom.push_back(tempInt);
      tempInt = 0;
      cout << endl;
      cout << "Please specify intermediate particle, if any ";
      cout << "(0 for none, else PDG code): ";
      cin >> tempInt;
      inter.push_back(tempInt);
      cout << endl;
      do {
        tempInt = 0;
        cout << "Please specify outgoing particle ";
        cout << "(jet=2212, else PDG code, exit with 99): ";
        cin >> tempInt;
        if (tempInt != 99) outgo.push_back(tempInt);
      } while(tempInt != 99);
      cout << endl;
    } else {
      cout << "LHE file not produced by SHERPA or MG/ME - ";
      cout << "Using default process and tms" << endl;
      incom.push_back(2212);
      incom.push_back(2212);
      inter.push_back(24);
      outgo.push_back(-11);
      outgo.push_back(12);
      tms = 10.;
    }
  }

  // Now store incoming, intermediate and outgoing
  // Set intermediate tags
  for(int i=0; i < int(inter.size()); ++i)
    hardIntermediate.push_back(inter[i]);

  // Set the incoming particle tags
  if (incom.size() != 2)
    cout << "Only two incoming particles allowed" << endl;
  else {
    hardIncoming1 = incom[0];
    hardIncoming2 = incom[1];
  }

  // Remember final state bosons
  int nBosons = 0;
  for(int i=0; i < int(outgo.size()); ++i)
    if ( (abs(outgo[i]) > 20 && abs(outgo[i]) <= 25) || outgo[i] == 2400)
      nBosons++;
  // Remember b-quark container
  int nBQuarks = 0;
  for(int i=0; i < int(outgo.size()); ++i)
    if ( outgo[i] == 5000)
      nBQuarks++;
  // Remember jet container
  int nJets = 0;
  for(int i=0; i < int(outgo.size()); ++i)
    if ( outgo[i] == 2212)
      nJets++;
  // Remember lepton container
  int nLeptons = 0;
  for(int i=0; i < int(outgo.size()); ++i)
    if ( outgo[i] == 1100)
      nLeptons++;
  // Remember lepton container
  int nNeutrinos = 0;
  for(int i=0; i < int(outgo.size()); ++i)
    if ( outgo[i] == 1200)
      nNeutrinos++;
  int nContainers = nLeptons + nNeutrinos + nJets + nBQuarks;

  // Set final particle identifiers
  if ( (outgo.size() - nBosons - nContainers)%2 == 1) {
    cout << "Only even number of outgoing particles allowed" << endl;
    for(int i=0; i < int(outgo.size()); ++i)
      cout << outgo[i] << endl;
  } else {

    // Push back particles / antiparticles
    for(int i=0; i < int(outgo.size()); ++i)
      if (outgo[i] > 0
        && outgo[i] != 2212
        && outgo[i] != 5000
        && outgo[i] != 1100
        && outgo[i] != 1200
        && outgo[i] != 2400
        && outgo[i] != 1000022)
        hardOutgoing2.push_back( outgo[i]);
      else if (outgo[i] < 0)
        hardOutgoing1.push_back( outgo[i]);

    // Save final state W-boson container as particle
    for(int i=0; i < int(outgo.size()); ++i)
      if ( outgo[i] == 2400)
        hardOutgoing2.push_back( outgo[i]);

    // Push back jets, distribute evenly amongst particles / antiparticles
    // Push back majorana particles, distribute evenly
    int iNow = 0;
    for(int i=0; i < int(outgo.size()); ++i)
      if ( (outgo[i] == 2212
        || outgo[i] == 5000
        || outgo[i] == 1200
        || outgo[i] == 1000022)
        && iNow%2 == 0 ){
        hardOutgoing2.push_back( outgo[i]);
        iNow++;
      } else if ( (outgo[i] == 2212
               || outgo[i] == 5000
               || outgo[i] == 1100
               || outgo[i] == 1000022)
               && iNow%2 == 1 ){
        hardOutgoing1.push_back( outgo[i]);
        iNow++;
      }
  }

  // Done
}

//--------------------------------------------------------------------------

// Function to translate a string specitying the core process into the
// internal notation
// Currently, the input string has to be in MadEvent notation

void HardProcess::translateProcessString( string process){

  vector <int> incom;
  vector <int> inter;
  vector <int> outgo;
  // Particle identifiers, ordered in such a way that e.g. the "u"
  // in a mu is not mistaken for an u quark
  int inParticleNumbers[] = {
        // Leptons
        -11,11,-12,12,-13,13,-14,14,-15,15,-16,16,
        // Jet container
        2212,2212,0,0,0,0,
        // Quarks
        -1,1,-2,2,-3,3,-4,4,-5,5,-6,6};
  string inParticleNamesMG[] =  {
        // Leptons
        "e+","e-","ve~","ve","mu+","mu-","vm~","vm","ta+","ta-","vt~","vt",
        // Jet container
        "p~","p","l+","l-","vl~","vl",
        // Quarks
        "d~","d","u~","u","s~","s","c~","c","b~","b","t~","t"};

  // Declare intermediate particle identifiers
  int interParticleNumbers[] = {
         // Electroweak gauge bosons
         22,23,-24,24,25,2400,
         // All squarks
        -1000001,1000001,-1000002,1000002,-1000003,1000003,-1000004,1000004,
        -1000005,1000005,-1000006,1000006,-2000001,2000001,-2000002,2000002,
        -2000003,2000003,-2000004,2000004,-2000005,2000005,-2000006,2000006,
         // Top quarks
         -6,6,
         // Dummy index as back-up
         0};
  // Declare names of intermediate particles
  string interParticleNamesMG[] = {
        // Electroweak gauge bosons
        "a","z","w-","w+","h","W",
         // All squarks
         "dl~","dl","ul~","ul","sl~","sl","cl~","cl","b1~","b1","t1~","t1",
         "dr~","dr","ur~","ur","sr~","sr","cr~","cr","b2~","b2","t2~","t2",
         // Top quarks
         "t~","t",
        // Dummy index as back-up
        "xx"};

  // Declare final state particle identifiers
  int outParticleNumbers[] = {
        // Leptons
        -11,11,-12,12,-13,13,-14,14,-15,15,-16,16,
        // Jet container and lepton containers
        2212,2212,0,0,0,0,1200,1100,5000,
        // All squarks
        -1000001,1000001,-1000002,1000002,-1000003,1000003,-1000004,1000004,
        -1000005,1000005,-1000006,1000006,-2000001,2000001,-2000002,2000002,
        -2000003,2000003,-2000004,2000004,-2000005,2000005,-2000006,2000006,
        // Quarks
        -1,1,-2,2,-3,3,-4,4,-5,5,-6,6,
        // SM uncoloured bosons
        22,23,-24,24,25,2400,
        // Neutralino in SUSY
        1000022};
  // Declare names of final state particles
  string outParticleNamesMG[] =  {
        // Leptons
        "e+","e-","ve~","ve","mu+","mu-","vm~","vm","ta+","ta-","vt~","vt",
        // Jet container and lepton containers
        "j~","j","l+","l-","vl~","vl","NEUTRINOS","LEPTONS","BQUARKS",
        // All squarks
        "dl~","dl","ul~","ul","sl~","sl","cl~","cl","b1~","b1","t1~","t1",
        "dr~","dr","ur~","ur","sr~","sr","cr~","cr","b2~","b2","t2~","t2",
        // Quarks
        "d~","d","u~","u","s~","s","c~","c","b~","b","t~","t",
        // SM uncoloured bosons
        "a","z","w-","w+","h","W",
        // Neutralino in SUSY
        "n1"};

  // Declare size of particle name arrays
  int nIn   = 30;
  int nInt  = 33;
  int nOut  = 64;

  // Start mapping user-defined particles onto particle ids.
  //string fullProc = "pp>{blaa,124}LEPTONS,NEUTRINOS";
  string fullProc = process;

  // Find user-defined hard process content
  // Count number of user particles
  int nUserParticles = 0;
  for(int n = fullProc.find("{", 0); n != int(string::npos);
          n = fullProc.find("{", n)) {
    nUserParticles++;
    n++;
  }
  // Cut user-defined particles from remaining process
  vector <string> userParticleStrings;
  for(int i =0; i < nUserParticles;++i) {
    int n = fullProc.find("{", 0);
    userParticleStrings.push_back(fullProc.substr(0,n));
    fullProc = fullProc.substr(n+1,fullProc.size());
  }
  // Cut remaining process string from rest
  if (nUserParticles > 0)
    userParticleStrings.push_back(
      fullProc.substr( 0, fullProc.find("}",0) ) );
  // Remove curly brackets and whitespace
  for(int i =0; i < int(userParticleStrings.size());++i) {
    while(userParticleStrings[i].find("{", 0) != string::npos)
      userParticleStrings[i].erase(userParticleStrings[i].begin()
                                  +userParticleStrings[i].find("{", 0));
    while(userParticleStrings[i].find("}", 0) != string::npos)
      userParticleStrings[i].erase(userParticleStrings[i].begin()
                                  +userParticleStrings[i].find("}", 0));
    while(userParticleStrings[i].find(" ", 0) != string::npos)
      userParticleStrings[i].erase(userParticleStrings[i].begin()
                                  +userParticleStrings[i].find(" ", 0));
  }
  // Convert particle numbers in user particle to integers
  vector<int>userParticleNumbers;
  if ( int(userParticleStrings.size()) > 1) {
    for( int i = 1; i < int(userParticleStrings.size()); ++i) {
      userParticleNumbers.push_back(
        atoi((char*)userParticleStrings[i].substr(
          userParticleStrings[i].find(",",0)+1,
          userParticleStrings[i].size()).c_str() ) );
    }
  }
  // Save remaining process string
  if (nUserParticles > 0)
    userParticleStrings.push_back(
      fullProc.substr(
        fullProc.find("}",0)+1, fullProc.size() ) );
  // Remove curly brackets and whitespace
  for( int i = 0; i < int(userParticleStrings.size()); ++i) {
    while(userParticleStrings[i].find("{", 0) != string::npos)
      userParticleStrings[i].erase(userParticleStrings[i].begin()
                                  +userParticleStrings[i].find("{", 0));
    while(userParticleStrings[i].find("}", 0) != string::npos)
      userParticleStrings[i].erase(userParticleStrings[i].begin()
                                  +userParticleStrings[i].find("}", 0));
    while(userParticleStrings[i].find(" ", 0) != string::npos)
      userParticleStrings[i].erase(userParticleStrings[i].begin()
                                  +userParticleStrings[i].find(" ", 0));
  }

  // Start mapping residual process string onto particle IDs.
  // Declare leftover process after user-defined particles have been converted
  string residualProc;
  if ( int(userParticleStrings.size()) > 1 )
    residualProc = userParticleStrings.front() + userParticleStrings.back();
  else
    residualProc = fullProc;

  // Remove comma separation
  while(residualProc.find(",", 0) != string::npos)
    residualProc.erase(residualProc.begin()+residualProc.find(",",0));

  // Count number of resonances
  int appearances = 0;
  for(int n = residualProc.find("(", 0); n != int(string::npos);
          n = residualProc.find("(", n)) {
    appearances++;
    n++;
  }

  // Cut string in incoming, resonance+decay and outgoing pieces
  vector <string> pieces;
  for(int i =0; i < appearances;++i) {
    int n = residualProc.find("(", 0);
    pieces.push_back(residualProc.substr(0,n));
    residualProc = residualProc.substr(n+1,residualProc.size());
  }
  // Cut last resonance from rest
  if (appearances > 0) {
    pieces.push_back( residualProc.substr(0,residualProc.find(")",0)) );
    pieces.push_back( residualProc.substr(
      residualProc.find(")",0)+1, residualProc.size()) );
  }

  // If the string was not cut into pieces, i.e. no resonance was
  // required, cut string using '>' as delimiter
  if (pieces.empty() ){
    appearances = 0;
    for(int n = residualProc.find(">", 0); n != int(string::npos);
            n = residualProc.find(">", n)) {
      appearances++;
      n++;
    }

    // Cut string in incoming and outgoing pieces
    for(int i =0; i < appearances;++i) {
      int n = residualProc.find(">", 0);
      pieces.push_back(residualProc.substr(0,n));
      residualProc = residualProc.substr(n+1,residualProc.size());
    }

    if (appearances == 1) pieces.push_back(residualProc);
    if (appearances > 1) {
      pieces.push_back( residualProc.substr(0,residualProc.find(">",0)) );
      pieces.push_back( residualProc.substr(
        residualProc.find(">",0)+1, residualProc.size()) );
    }
  }

  // Get incoming particles
  for(int i=0; i < nIn; ++i) {
    for(int n = pieces[0].find(inParticleNamesMG[i], 0);
           n != int(string::npos);
           n = pieces[0].find(inParticleNamesMG[i], n)) {
      incom.push_back(inParticleNumbers[i]);
      pieces[0].erase(pieces[0].begin()+n,
                      pieces[0].begin()+n+inParticleNamesMG[i].size());
      n=0;
    }
  }

  // Check intermediate resonances and decay products
  for(int i =1; i < int(pieces.size()); ++i){
    // Seperate strings into intermediate and outgoing, if not already done
    int k = pieces[i].find(">", 0);

    string intermediate = (pieces[i].find(">", 0) != string::npos) ?
                           pieces[i].substr(0,k) : "";
    string outgoing = (pieces[i].find(">", 0) != string::npos) ?
                       pieces[i].substr(k+1,pieces[i].size()) : pieces[i];

    // Get intermediate particles
    for(int j=0; j < nInt; ++j) {
      for(int n = intermediate.find(interParticleNamesMG[j], 0);
             n != int(string::npos);
             n = intermediate.find(interParticleNamesMG[j], n)) {
        inter.push_back(interParticleNumbers[j]);
        intermediate.erase(intermediate.begin()+n,
                    intermediate.begin()+n+interParticleNamesMG[j].size());
        n=0;
      }
    }

    // Get outgoing particles
    for(int j=0; j < nOut; ++j) {
      for(int n = outgoing.find(outParticleNamesMG[j], 0);
             n != int(string::npos);
             n = outgoing.find(outParticleNamesMG[j], n)) {
        outgo.push_back(outParticleNumbers[j]);
        outgoing.erase(outgoing.begin()+n,
                       outgoing.begin()+n+outParticleNamesMG[j].size());
        n=0;
      }
    }

    // For arbitrary or non-existing intermediate, remember zero for each
    // two outgoing particles, without bosons.
    if (inter.empty()) {

      // For final state bosons, bookkeep the final state boson as
      // intermediate as well.
      int nBosons = 0;
      for(int l=0; l < int(outgo.size()); ++l)
        if ( (abs(outgo[l]) > 20 && abs(outgo[l]) <= 25) || outgo[l] == 2400)
          nBosons++;

      int nZeros = (outgo.size() - nBosons)/2;
      for(int l=0; l < nZeros; ++l)
        inter.push_back(0);
    }

    // For final state bosons, bookkeep the final state boson as
    // intermediate as well.
    for(int l=0; l < int(outgo.size()); ++l)
      if ( (abs(outgo[l]) > 20 && abs(outgo[l]) <= 25) || outgo[l] == 2400)
        inter.push_back(outgo[l]);

  }

  // Now store incoming, intermediate and outgoing
  // Set intermediate tags
  for(int i=0; i < int(inter.size()); ++i)
    hardIntermediate.push_back(inter[i]);

  // Set the incoming particle tags
  if (incom.size() != 2)
    cout << "Only two incoming particles allowed" << endl;
  else {
    hardIncoming1 = incom[0];
    hardIncoming2 = incom[1];
  }

  // Now store final particle identifiers
  // Start with user-defined particles.
  for( int i = 0; i < int(userParticleNumbers.size()); ++i)
    if (userParticleNumbers[i] > 0) {
      hardOutgoing2.push_back( userParticleNumbers[i]);
      hardIntermediate.push_back(0);
      // For non-existing intermediate, remember zero.
    } else if (userParticleNumbers[i] < 0) {
      hardOutgoing1.push_back( userParticleNumbers[i]);
      // For non-existing intermediate, remember zero.
      hardIntermediate.push_back(0);
    }

  // Push back particles / antiparticles
  for(int i=0; i < int(outgo.size()); ++i)
    if (outgo[i] > 0
      && outgo[i] != 2212
      && outgo[i] != 5000
      && outgo[i] != 1100
      && outgo[i] != 1200
      && outgo[i] != 2400
      && outgo[i] != 1000022)
      hardOutgoing2.push_back( outgo[i]);
    else if (outgo[i] < 0)
      hardOutgoing1.push_back( outgo[i]);

  // Save final state W-boson container as particle
  for(int i=0; i < int(outgo.size()); ++i)
    if ( outgo[i] == 2400)
      hardOutgoing2.push_back( outgo[i]);

  // Push back jets, distribute evenly among particles / antiparticles
  // Push back majorana particles, distribute evenly
  int iNow = 0;
  for(int i=0; i < int(outgo.size()); ++i)
    if ( (outgo[i] == 2212
      || outgo[i] == 5000
      || outgo[i] == 1200
      || outgo[i] == 1000022)
      && iNow%2 == 0 ){
      hardOutgoing2.push_back( outgo[i]);
      iNow++;
    } else if ( (outgo[i] == 2212
             || outgo[i] == 5000
             || outgo[i] == 1100
             || outgo[i] == 1000022)
             && iNow%2 == 1 ){
      hardOutgoing1.push_back( outgo[i]);
      iNow++;
    }

  // Done
}

//--------------------------------------------------------------------------

// Function to check if the candidates stored in Pos1 and Pos2, together with
// a proposed candidate iPos are allowed.

bool HardProcess::allowCandidates(int iPos, vector<int> Pos1,
  vector<int> Pos2, const Event& event){

  bool allowed = true;

  // Find colour-partner of new candidate
  int type = (event[iPos].col() > 0) ? 1 : (event[iPos].acol() > 0) ? -1 : 0;

  if (type == 0) return true;

  if (type == 1){
    int col = event[iPos].col();
    int iPartner = 0;
    for(int i=0; i < int(event.size()); ++i)
      if ( i != iPos
        && (( event[i].isFinal() && event[i].acol() == col)
          ||( event[i].status() == -21 && event[i].col() == col) ))
      iPartner = i;

    vector<int> partners;
    for(int i=0; i < int(event.size()); ++i)
      for(int j=0; j < int(Pos1.size()); ++j)
        if ( Pos1[j] != 0 && i != Pos1[j] && event[Pos1[j]].colType() != 0
        && (( event[i].isFinal() && event[i].col() == event[Pos1[j]].acol())
          ||( event[i].status() == -21
           && event[i].acol() == event[Pos1[j]].acol()) ))
         partners.push_back(i);

    // Never allow equal initial partners!
    if (event[iPartner].status() == -21){
      for(int i=0; i < int(partners.size()); ++i)
        if ( partners[i] == iPartner)
          allowed = false;
    }

  } else {
    int col = event[iPos].acol();
    int iPartner = 0;
    for(int i=0; i < int(event.size()); ++i)
      if ( i != iPos
        && (( event[i].isFinal() && event[i].col()  == col)
          ||(!event[i].isFinal() && event[i].acol() == col) ))
      iPartner = i;

    vector<int> partners;
    for(int i=0; i < int(event.size()); ++i)
      for(int j=0; j < int(Pos2.size()); ++j)
        if ( Pos2[j] != 0 && i != Pos2[j] && event[Pos2[j]].colType() != 0
        && (( event[i].isFinal() && event[i].acol() == event[Pos2[j]].col())
          ||( event[i].status() == -21
           && event[i].col() == event[Pos2[j]].col()) ))
         partners.push_back(i);

    // Never allow equal initial partners!
    if (event[iPartner].status() == -21){
      for(int i=0; i < int(partners.size()); ++i){
        if ( partners[i] == iPartner)
          allowed = false;
      }
    }

  }

  return allowed;

}

//--------------------------------------------------------------------------

// Function to identify the hard subprocess in the current event

void HardProcess::storeCandidates( const Event& event, string process){

  // Store the reference event
  state.clear();
  state = event;

  // Local copy of intermediate bosons
  vector<int> intermediates;
  for(int i =0; i < int(hardIntermediate.size());++i)
    intermediates.push_back( hardIntermediate[i]);

  // Local copy of outpoing partons
  vector<int> outgoing1;
  for(int i =0; i < int(hardOutgoing1.size());++i)
    outgoing1.push_back( hardOutgoing1[i]);
  vector<int> outgoing2;
  for(int i =0; i < int(hardOutgoing2.size());++i)
    outgoing2.push_back( hardOutgoing2[i]);

  // Clear positions of intermediate and outgoing particles
  PosIntermediate.resize(0);
  PosOutgoing1.resize(0);
  PosOutgoing2.resize(0);
  for(int i =0; i < int(hardIntermediate.size());++i)
    PosIntermediate.push_back(0);
  for(int i =0; i < int(hardOutgoing1.size());++i)
    PosOutgoing1.push_back(0);
  for(int i =0; i < int(hardOutgoing2.size());++i)
    PosOutgoing2.push_back(0);

  // For QCD dijet or e+e- > jets hard process, do not store any candidates,
  // as to not discrimintate clusterings
  if (  process.compare("pp>jj") == 0
    || process.compare("e+e->jj") == 0
    || process.compare("e+e->(z>jj)") == 0 ){
    for(int i =0; i < int(hardOutgoing1.size());++i)
      PosOutgoing1[i] = 0;
    for(int i =0; i < int(hardOutgoing2.size());++i)
      PosOutgoing2[i] = 0;
    // Done
    return;
  }

  // Initialise vector of particles that were already identified as
  // hard process particles
  vector<int> iPosChecked;

  // If the hard process is specified only by containers, then add all
  // particles matching with the containers to the hard process.
  bool hasOnlyContainers = true;
  for(int i =0; i < int(hardOutgoing1.size());++i)
    if (  hardOutgoing1[i] != 1100
      && hardOutgoing1[i] != 1200
      && hardOutgoing1[i] != 5000)
      hasOnlyContainers = false;
  for(int i =0; i < int(hardOutgoing2.size());++i)
    if (  hardOutgoing2[i] != 1100
      && hardOutgoing2[i] != 1200
      && hardOutgoing2[i] != 5000)
      hasOnlyContainers = false;

  if (hasOnlyContainers){

    PosOutgoing1.resize(0);
    PosOutgoing2.resize(0);

    // Try to find all unmatched hard process leptons.
    // Loop through event to find outgoing lepton
    for(int i=0; i < int(event.size()); ++i){

      // Skip non-final particles
      if ( !event[i].isFinal() ) continue;

      // Skip all particles that have already been identified
      bool skip = false;
      for(int k=0; k < int(iPosChecked.size()); ++k){
        if (i == iPosChecked[k])
          skip = true;
      }
      if (skip) continue;

      for(int j=0; j < int(outgoing2.size()); ++j){

        // If the particle matches an outgoing neutrino, save it
        if ( outgoing2[j] == 1100
          && ( event[i].idAbs() == 11
            || event[i].idAbs() == 13
            || event[i].idAbs() == 15) ){
          PosOutgoing2.push_back(i);
          iPosChecked.push_back(i);
        }

        // If the particle matches an outgoing lepton, save it
        if ( outgoing2[j] == 1200
          && ( event[i].idAbs() == 12
            || event[i].idAbs() == 14
            || event[i].idAbs() == 16) ){
          PosOutgoing2.push_back(i);
          iPosChecked.push_back(i);
        }

        // If the particle matches an outgoing b-quark, save it
        if ( outgoing2[j] == 5000 && event[i].idAbs() == 5 ){
          PosOutgoing2.push_back(i);
          iPosChecked.push_back(i);
        }

      }

      // Skip all particles that have already been identified
      skip = false;
      for(int k=0; k < int(iPosChecked.size()); ++k){
        if (i == iPosChecked[k])
          skip = true;
      }
      if (skip) continue;

      for(int j=0; j < int(outgoing1.size()); ++j){

        // If the particle matches an outgoing neutrino, save it
        if ( outgoing1[j] == 1100
          && ( event[i].idAbs() == 11
            || event[i].idAbs() == 13
            || event[i].idAbs() == 15) ){
          PosOutgoing1.push_back(i);
          iPosChecked.push_back(i);
        }

        // If the particle matches an outgoing lepton, save it
        if ( outgoing1[j] == 1200
          && ( event[i].idAbs() == 12
            || event[i].idAbs() == 14
            || event[i].idAbs() == 16) ){
          PosOutgoing1.push_back(i);
          iPosChecked.push_back(i);
        }

        // If the particle matches an outgoing b-quark, save it
        if ( outgoing1[j] == 5000 && event[i].idAbs() == 5 ){
          PosOutgoing1.push_back(i);
          iPosChecked.push_back(i);
        }

      }
    }

    // Done
    return;
  }

  // Now begin finding candidates when not only containers are used.

  // First try to find final state bosons
  for(int i=0; i < int(intermediates.size()); ++i){

    // Do nothing if the intermediate boson is absent
    if (intermediates[i] == 0) continue;

    // Do nothing if this boson does not match any final state boson
    bool matchesFinalBoson = false;
    for(int j =0; j< int(outgoing1.size()); ++j){
      if ( intermediates[i] == outgoing1[j] )
        matchesFinalBoson = true;
    }
    for(int j =0; j< int(outgoing2.size()); ++j){
      if ( intermediates[i] == outgoing2[j] )
        matchesFinalBoson = true;
    }
    if (!matchesFinalBoson) continue;

    // Loop through event
    for(int j=0; j < int(event.size()); ++j) {

      // Skip all particles that have already been identified
      bool skip = false;
      for(int m=0; m < int(iPosChecked.size()); ++m)
        if (j == iPosChecked[m]) skip = true;
      if (skip) continue;

      // If the particle has a requested intermediate id, check if
      // if is a final state boson
      if ( (event[j].id() == intermediates[i])
        ||(event[j].idAbs() == 24 && intermediates[i] == 2400) ) {

        PosIntermediate[i] = j;
        intermediates[i] = 0;
        // Be careful only to replace one index at a time!
        bool indexSet = false;

        for(int k=0; k < int(outgoing1.size()); ++k) {
          if (event[j].id() == outgoing1[k] && !indexSet){
            PosOutgoing1[k] = j;
            outgoing1[k] = 99;
            indexSet = true;
          }
        }

        for(int k=0; k < int(outgoing2.size()); ++k) {
          if (event[j].id() == outgoing2[k] && !indexSet){
            PosOutgoing2[k] = j;
            outgoing2[k] = 99;
            indexSet = true;
          }
        }

        // Check for W-boson container
        for(int k=0; k < int(outgoing2.size()); ++k) {
          if (event[j].idAbs() == 24 && outgoing2[k] == 2400 && !indexSet ){
            PosOutgoing2[k] = j;
            outgoing2[k] = 99;
            indexSet = true;
          }
        }

        iPosChecked.push_back(j);

      }
    }
  }

  // Second try to find particles coupled to intermediate bosons
  for(int i=0; i < int(intermediates.size()); ++i){

    // Do nothing if the intermediate boson is absent
    if (intermediates[i] == 0) continue;

    // Loop through event
    for(int j=0; j < int(event.size()); ++j) {
      // If the particle has a requested intermediate id, check if
      // daughters are hard process particles
      if ( (event[j].id() == intermediates[i])
        ||(event[j].idAbs() == 24 && intermediates[i] == 2400) ) {
        // If this particle is a potential intermediate
        PosIntermediate[i] = j;
        intermediates[i] = 0;
        // If id's of daughters are good, store position
        int iPos1 = event[j].daughter1();
        int iPos2 = event[j].daughter2();

        // Loop through daughters to check if these contain some hard
        // outgoing particles
        for( int k=iPos1; k <= iPos2; ++k){
          int id = event[k].id();

          // Skip all particles that have already been identified
          bool skip = false;
          for(int m=0; m < int(iPosChecked.size()); ++m)
            if (k == iPosChecked[m]) skip = true;
          if (skip) continue;

          // Check if daughter is hard outgoing particle
          for(int l=0; l < int(outgoing2.size()); ++l)
            if ( outgoing2[l] != 99 ){
                // Found particle id
              if (id == outgoing2[l]
                // Found jet
                || (id > 0 && abs(id) < 10 && outgoing2[l] == 2212) ){
                // Store position
                PosOutgoing2[l] = k;
                // Remove the matched particle from the list
                outgoing2[l] = 99;
                iPosChecked.push_back(k);
                break;
              }

            }

          // Check if daughter is hard outgoing antiparticle
          for(int l=0; l < int(outgoing1.size()); ++l)
            if ( outgoing1[l] != 99 ){
                // Found particle id
              if (id == outgoing1[l]
                // Found jet
                || (id < 0 && abs(id) < 10 && outgoing1[l] == 2212) ){
                // Store position
                PosOutgoing1[l] = k;
                // Remove the matched particle from the list
                outgoing1[l] = 99;
                iPosChecked.push_back(k);
                break;
            }

          }

        } // End loop through daughters
      } // End if ids match
    } // End loop through event
  } // End loop though requested intermediates

  // If all outgoing particles were found, done
  bool done = true;
  for(int i=0; i < int(outgoing1.size()); ++i)
    if (outgoing1[i] != 99)
      done = false;
  for(int i=0; i < int(outgoing2.size()); ++i)
    if (outgoing2[i] != 99)
      done = false;
  // Return
  if (done) return;

  // Leptons not associated with resonance are allowed.
  // Try to find all unmatched hard process leptons.
  // Loop through event to find outgoing lepton
  for(int i=0; i < int(event.size()); ++i){
    // Skip non-final particles and final partons
    if ( !event[i].isFinal() || event[i].colType() != 0)
      continue;
    // Skip all particles that have already been identified
    bool skip = false;
    for(int k=0; k < int(iPosChecked.size()); ++k){
      if (i == iPosChecked[k])
        skip = true;
    }
    if (skip) continue;

    // Check if any hard outgoing leptons remain
    for(int j=0; j < int(outgoing2.size()); ++j){
      // Do nothing if this particle has already be found,
      // or if this particle is a jet or quark
      if (  outgoing2[j] == 99
        || outgoing2[j] == 2212
        || abs(outgoing2[j]) < 10)
        continue;

      // If the particle matches an outgoing lepton, save it
      if (  event[i].id() == outgoing2[j] ){
        PosOutgoing2[j] = i;
        outgoing2[j] = 99;
        iPosChecked.push_back(i);
      }
    }

    // Check if any hard outgoing antileptons remain
    for(int j=0; j < int(outgoing1.size()); ++j){
      // Do nothing if this particle has already be found,
      // or if this particle is a jet or quark
      if (  outgoing1[j] == 99
        || outgoing1[j] == 2212
        || abs(outgoing1[j]) < 10)
        continue;

      // If the particle matches an outgoing lepton, save it
      if (event[i].id() == outgoing1[j] ){
        PosOutgoing1[j] = i;
        outgoing1[j] = 99;
        iPosChecked.push_back(i);
      }
    }
  }

  multimap<int,int> out2copy;
  for(int i=0; i < int(event.size()); ++i)
    for(int j=0; j < int(outgoing2.size()); ++j)
      // Do nothing if this particle has already be found,
      // or if this particle is a jet.
      if ( outgoing2[j] != 99
        && outgoing2[j] != 2212
        && ( abs(outgoing2[j]) < 10
          || (abs(outgoing2[j]) > 1000000 && abs(outgoing2[j]) < 1000010)
          || (abs(outgoing2[j]) > 2000000 && abs(outgoing2[j]) < 2000010)
          || abs(outgoing2[j]) == 1000021 )
        && event[i].isFinal()
        && event[i].id() == outgoing2[j] ){
        out2copy.insert(make_pair(j, i));
      }

  multimap<int,int> out1copy;
  for(int i=0; i < int(event.size()); ++i)
    for(int j=0; j < int(outgoing1.size()); ++j)
      // Do nothing if this particle has already be found,
      // or if this particle is a jet.
      if ( outgoing1[j] != 99
        && outgoing1[j] != 2212
        && ( abs(outgoing1[j]) < 10
          || (abs(outgoing1[j]) > 1000000 && abs(outgoing1[j]) < 1000010)
          || (abs(outgoing1[j]) > 2000000 && abs(outgoing1[j]) < 2000010)
          || abs(outgoing2[j]) == 1000021 )
        && event[i].isFinal()
        && event[i].id() == outgoing1[j] ){
        out1copy.insert(make_pair(j, i));
      }

  if ( out1copy.size() >  out2copy.size()){

    // In case the index of the multimap is filled twice, make sure not to
    // arbitrarily overwrite set values.
    vector<int> indexWasSet;
    for ( multimap<int, int>::iterator it = out2copy.begin();
      it != out2copy.end(); ++it ) {
      if ( allowCandidates(it->second, PosOutgoing1, PosOutgoing2, event) ){

        // Skip all particles that have already been identified
        bool skip = false;
        for(int k=0; k < int(iPosChecked.size()); ++k)
          if (it->second == iPosChecked[k]) skip = true;
        // Skip all indices that have already been identified
        for(int k=0; k < int(indexWasSet.size()); ++k)
          if (it->first == indexWasSet[k]) skip = true;
        if (skip) continue;

        // Save parton
        PosOutgoing2[it->first] = it->second;
        // remove entry form lists
        outgoing2[it->first] = 99;
        iPosChecked.push_back(it->second);
        indexWasSet.push_back(it->first);
      }
    }

    indexWasSet.resize(0);
    for ( multimap<int, int>::iterator it = out1copy.begin();
      it != out1copy.end(); ++it ) {
      if ( allowCandidates(it->second, PosOutgoing1, PosOutgoing2, event) ){

        // Skip all particles that have already been identified
        bool skip = false;
        for(int k=0; k < int(iPosChecked.size()); ++k)
          if (it->second == iPosChecked[k]) skip = true;
        // Skip all indices that have already been identified
        for(int k=0; k < int(indexWasSet.size()); ++k)
          if (it->first == indexWasSet[k]) skip = true;
        if (skip) continue;

        // Save parton
        PosOutgoing1[it->first] = it->second;
        // remove entry form lists
        outgoing1[it->first] = 99;
        iPosChecked.push_back(it->second);
        indexWasSet.push_back(it->first);
      }
    }

  } else {

    // In case the index of the multimap is filled twice, make sure not to
    // arbitraryly overwrite set values.
    vector<int> indexWasSet;
    for ( multimap<int, int>::iterator it = out1copy.begin();
      it != out1copy.end(); ++it ) {
      if ( allowCandidates(it->second, PosOutgoing1, PosOutgoing2, event) ){

        // Skip all particles that have already been identified
        bool skip = false;
        for(int k=0; k < int(iPosChecked.size()); ++k)
          if (it->second == iPosChecked[k]) skip = true;
        // Skip all indices that have already been identified
        for(int k=0; k < int(indexWasSet.size()); ++k)
          if (it->first == indexWasSet[k]) skip = true;
        if (skip) continue;

        // Save parton
        PosOutgoing1[it->first] = it->second;
        // remove entry form lists
        outgoing1[it->first] = 99;
        iPosChecked.push_back(it->second);
        indexWasSet.push_back(it->first);
      }
    }

    indexWasSet.resize(0);
    for ( multimap<int, int>::iterator it = out2copy.begin();
      it != out2copy.end(); ++it ) {
      if ( allowCandidates(it->second, PosOutgoing1, PosOutgoing2, event) ){

        // Skip all particles that have already been identified
        bool skip = false;
        for(int k=0; k < int(iPosChecked.size()); ++k)
          if (it->second == iPosChecked[k]) skip = true;
        // Skip all indices that have already been identified
        for(int k=0; k < int(indexWasSet.size()); ++k)
          if (it->first == indexWasSet[k]) skip = true;
        if (skip) continue;

        // Save parton
        PosOutgoing2[it->first] = it->second;
        // remove entry form lists
        outgoing2[it->first] = 99;
        iPosChecked.push_back(it->second);
        indexWasSet.push_back(it->first);
      }
    }
  }

  // It sometimes happens that MadEvent does not put a
  // heavy coloured resonance into the LHE file, even if requested.
  // This means that the decay products of this resonance need to be
  // found separately.
  // Loop through event to find hard process (anti)quarks
  for(int i=0; i < int(event.size()); ++i){

    // Skip non-final particles and final partons
    if ( !event[i].isFinal() || event[i].colType() == 0)
      continue;

    // Skip all particles that have already been identified
    bool skip = false;
    for(int k=0; k < int(iPosChecked.size()); ++k){
      if (i == iPosChecked[k])
        skip = true;
    }
    if (skip) continue;

    // Check if any hard outgoing quarks remain
    for(int j=0; j < int(outgoing2.size()); ++j){
      // Do nothing if this particle has already be found,
      // or if this particle is a jet, lepton container or lepton
      if (  outgoing2[j] == 99
        || outgoing2[j] == 2212
        || abs(outgoing2[j]) > 10)
        continue;
      // If the particle matches an outgoing quark, save it
      if (event[i].id() == outgoing2[j]){
        // Save parton
        PosOutgoing2[j] = i;
        // remove entry form lists
        outgoing2[j] = 99;
        iPosChecked.push_back(i);
      }
    }

    // Check if any hard outgoing antiquarks remain
    for(int j=0; j < int(outgoing1.size()); ++j){
      // Do nothing if this particle has already be found,
      // or if this particle is a jet, lepton container or lepton
      if (  outgoing1[j] == 99
        || outgoing1[j] == 2212
        || abs(outgoing1[j]) > 10)
        continue;
      // If the particle matches an outgoing antiquark, save it
      if (event[i].id() == outgoing1[j]){
        // Save parton
        PosOutgoing1[j] = i;
        // Remove parton from list
        outgoing1[j] = 99;
        iPosChecked.push_back(i);
      }
    }
  }

  // Done
}

//--------------------------------------------------------------------------

// Function to check if the particle event[iPos] matches any of
// the stored outgoing particles of the hard subprocess

bool HardProcess::matchesAnyOutgoing(int iPos, const Event& event){

  // Match quantum numbers of any first outgoing candidate
  bool matchQN1 = false;
  // Match quantum numbers of any second outgoing candidate
  bool matchQN2 = false;
  // Match parton in the hard process,
  // or parton from decay of electroweak boson in hard process,
  // or parton from decay of electroweak boson from decay of top
  bool matchHP = false;

  // Check outgoing candidates
  for(int i=0; i < int(PosOutgoing1.size()); ++i)
    // Compare particle properties
    if ( event[iPos].id()         == state[PosOutgoing1[i]].id()
     && event[iPos].colType()    == state[PosOutgoing1[i]].colType()
     && event[iPos].chargeType() == state[PosOutgoing1[i]].chargeType()
     && ( ( event[iPos].col() > 0
         && event[iPos].col() == state[PosOutgoing1[i]].col())
       || ( event[iPos].acol() > 0
         && event[iPos].acol() == state[PosOutgoing1[i]].acol()))
     && event[iPos].charge()     == state[PosOutgoing1[i]].charge() )
      matchQN1 = true;

  // Check outgoing candidates
  for(int i=0; i < int(PosOutgoing2.size()); ++i)
    // Compare particle properties
    if ( event[iPos].id()         == state[PosOutgoing2[i]].id()
     && event[iPos].colType()    == state[PosOutgoing2[i]].colType()
     && event[iPos].chargeType() == state[PosOutgoing2[i]].chargeType()
     && ( ( event[iPos].col() > 0
         && event[iPos].col() == state[PosOutgoing2[i]].col())
       || ( event[iPos].acol() > 0
         && event[iPos].acol() == state[PosOutgoing2[i]].acol()))
     && event[iPos].charge()     == state[PosOutgoing2[i]].charge() )
      matchQN2 = true;

  // Check if maps to hard process:
  // Check that particle is in hard process
  if ( event[iPos].mother1()*event[iPos].mother2() == 12
      // Or particle has taken recoil from first splitting
      || (  event[iPos].status() == 44
         && event[event[iPos].mother1()].mother1()
           *event[event[iPos].mother1()].mother2() == 12 )
      // Or particle has on-shell resonace as mother
      || (  event[iPos].status() == 23
         && event[event[iPos].mother1()].mother1()
           *event[event[iPos].mother1()].mother2() == 12 )
      // Or particle has on-shell resonace as mother,
      // which again has and on-shell resonance as mother
      || (  event[iPos].status() == 23
         && event[event[iPos].mother1()].status() == -22
         && event[event[event[iPos].mother1()].mother1()].status() == -22
         && event[event[event[iPos].mother1()].mother1()].mother1()
           *event[event[event[iPos].mother1()].mother1()].mother2() == 12 ) )
      matchHP = true;

  // Done
  return ( matchHP && (matchQN1 || matchQN2) );

}


//--------------------------------------------------------------------------

// Function to check if instead of the particle event[iCandidate], another
// particle could serve as part of the hard process. Assumes that iCandidate
// is already stored as part of the hard process.

bool HardProcess::findOtherCandidates(int iPos, const Event& event,
    bool doReplace){

  // Return value
  bool foundCopy = false;

  // Save stored candidates' properties.
  int id  = event[iPos].id();
  int col = event[iPos].col();
  int acl = event[iPos].acol();

  // If the particle's mother is an identified intermediate resonance,
  // then do not attempt any replacement.
  int iMoth1 = event[iPos].mother1();
  int iMoth2 = event[iPos].mother2();
  if ( iMoth1 > 0 && iMoth2 == 0 ) {
    bool hasIdentifiedMother = false;
    for(int i=0; i < int(PosIntermediate.size()); ++i)
      // Compare particle properties
      if ( event[iMoth1].id()         == state[PosIntermediate[i]].id()
        && event[iMoth1].colType()    == state[PosIntermediate[i]].colType()
        && event[iMoth1].chargeType() == state[PosIntermediate[i]].chargeType()
        && event[iMoth1].col()        == state[PosIntermediate[i]].col()
        && event[iMoth1].acol()       == state[PosIntermediate[i]].acol()
        && event[iMoth1].charge()     == state[PosIntermediate[i]].charge() )
         hasIdentifiedMother = true;
    if(hasIdentifiedMother && event[iMoth1].id() != id) return false;
  }

  // Find candidate amongst the already stored ME process candidates.
  vector<int> candidates1;
  vector<int> candidates2;
  // Check outgoing candidates
  for(int i=0; i < int(PosOutgoing1.size()); ++i)
    // Compare particle properties
    if ( id  == state[PosOutgoing1[i]].id()
      && col == state[PosOutgoing1[i]].col()
      && acl == state[PosOutgoing1[i]].acol() )
      candidates1.push_back(i);
  // Check outgoing candidates
  for(int i=0; i < int(PosOutgoing2.size()); ++i)
    // Compare particle properties
    if ( id  == state[PosOutgoing2[i]].id()
      && col == state[PosOutgoing2[i]].col()
      && acl == state[PosOutgoing2[i]].acol() )
      candidates2.push_back(i);

  // If more / less than one stored candidate for iPos has been found, exit.
  if ( candidates1.size() + candidates2.size() != 1 ) return false;

  // Now check for other allowed candidates.
  map<int,int> further1;
  for(int i=0; i < int(state.size()); ++i)
    for(int j=0; j < int(PosOutgoing1.size()); ++j)
      // Do nothing if this particle has already be found,
      // or if this particle is a jet, lepton container or lepton
      if ( state[i].isFinal()
        && i != PosOutgoing1[j]
        && state[PosOutgoing1[j]].id() == id
        && state[i].id() == id ){
        // Declare vector of already existing candiates.
        vector<int> newPosOutgoing1;
        for(int k=0; k < int(PosOutgoing1.size()); ++k)
          if ( k != j ) newPosOutgoing1.push_back( PosOutgoing1[k] );
        // If allowed, remember replacment parton.
        if ( allowCandidates(i, newPosOutgoing1, PosOutgoing2, state) )
          further1.insert(make_pair(j, i));
      }

  // Now check for other allowed candidates.
  map<int,int> further2;
  for(int i=0; i < int(state.size()); ++i)
    for(int j=0; j < int(PosOutgoing2.size()); ++j)
      // Do nothing if this particle has already be found,
      // or if this particle is a jet, lepton container or lepton
      if ( state[i].isFinal()
        && i != PosOutgoing2[j]
        && state[PosOutgoing2[j]].id() == id
        && state[i].id() == id ){
        // Declare vector of already existing candidates.
        vector<int> newPosOutgoing2;
        for(int k=0; k < int(PosOutgoing2.size()); ++k)
          if ( k != j ) newPosOutgoing2.push_back( PosOutgoing2[k] );
        // If allowed, remember replacment parton.
        if ( allowCandidates(i, PosOutgoing1, newPosOutgoing2, state) )
          further2.insert(make_pair(j, i));
      }

  // Remove all hard process particles that would be counted twice.
  map<int,int>::iterator it2 = further2.begin();
  while(it2 != further2.end()) {
    bool remove = false;
    for(int j=0; j < int(PosOutgoing2.size()); ++j)
      if (it2->second == PosOutgoing2[j] ) remove = true;
    if ( remove ) further2.erase(it2++);
    else ++it2;
  }
  map<int,int>::iterator it1 = further1.begin();
  while(it1 != further1.end()) {
    bool remove = false;
    for(int j=0; j < int(PosOutgoing1.size()); ++j)
      if (it1->second == PosOutgoing1[j] ) remove = true;
    if ( remove ) further1.erase(it1++);
    else ++it1;
  }

  // Decide of a replacment candidate has been found.
  foundCopy = (doReplace)
            ? exchangeCandidates(candidates1, candidates2, further1, further2)
            : (further1.size() + further2.size() > 0);

  // Done
  return foundCopy;

}

//--------------------------------------------------------------------------

// Function to exchange hard process candidates.

bool HardProcess::exchangeCandidates( vector<int> candidates1,
    vector<int> candidates2, map<int,int> further1, map<int,int> further2) {

  int nOld1 = candidates1.size();
  int nOld2 = candidates2.size();
  int nNew1 = further1.size();
  int nNew2 = further2.size();
  bool exchanged = false;
  // Replace, if one-to-one correspondence exists.
  if ( nOld1 == 1 && nOld2 == 0 && nNew1 == 1 && nNew2 == 0){
    PosOutgoing1[further1.begin()->first] = further1.begin()->second;
    exchanged = true;
  } else if ( nOld1 == 0 && nOld2 == 1 && nNew1 == 0 && nNew2 == 1){
    PosOutgoing2[further2.begin()->first] = further2.begin()->second;
    exchanged = true;
  // Else simply swap with the first candidate.
  } else if ( nNew1 >  1 && nNew2 == 0 ) {
    PosOutgoing1[further1.begin()->first] = further1.begin()->second;
    exchanged = true;
  } else if ( nNew1 == 0 && nNew2 >  0 ) {
    PosOutgoing2[further2.begin()->first] = further2.begin()->second;
    exchanged = true;
  }

  // Done
  return exchanged;

}

//--------------------------------------------------------------------------

// Function to get the number of coloured final state partons in the
// hard process

int HardProcess::nQuarksOut(){
  int nFin =0;
  for(int i =0; i< int(hardOutgoing1.size()); ++i){
    if (hardOutgoing1[i] == 2212 || abs(hardOutgoing1[i]) < 10) nFin++;
  }
  for(int i =0; i< int(hardOutgoing2.size()); ++i){
    if (hardOutgoing2[i] == 2212 || abs(hardOutgoing2[i]) < 10) nFin++;
  }
  // For very loose hard process definition, check number of hard process
  // b-quarks explicitly.
  for(int i =0; i< int(hardOutgoing1.size()); ++i)
    if (hardOutgoing1[i] == 5000)
      for(int j =0; j< int(PosOutgoing1.size()); ++j)
        if (state[PosOutgoing1[j]].idAbs() == 5)
          nFin++;
  for(int i =0; i< int(hardOutgoing2.size()); ++i)
    if (hardOutgoing2[i] == 5000)
      for(int j =0; j< int(PosOutgoing2.size()); ++j)
        if (state[PosOutgoing2[j]].idAbs() == 5)
          nFin++;
  return nFin;
}

//--------------------------------------------------------------------------

// Function to get the number of uncoloured final state particles in the
// hard process

int HardProcess::nLeptonOut(){
  int nFin =0;
  for(int i =0; i< int(hardOutgoing1.size()); ++i){
    if (abs(hardOutgoing1[i]) > 10 && abs(hardOutgoing1[i]) < 20) nFin++;
    // Bookkeep MSSM neutralinos as leptons
    if (abs(hardOutgoing1[i]) == 1000022) nFin++;
    // Bookkeep sleptons as leptons
    if ( abs(hardOutgoing1[i]) == 1000011 || abs(hardOutgoing1[i]) == 2000011
      || abs(hardOutgoing1[i]) == 1000013 || abs(hardOutgoing1[i]) == 2000013
      || abs(hardOutgoing1[i]) == 1000015 || abs(hardOutgoing1[i]) == 2000015)
      nFin++;
  }
  for(int i =0; i< int(hardOutgoing2.size()); ++i){
    if (abs(hardOutgoing2[i]) > 10 && abs(hardOutgoing2[i]) < 20) nFin++;
    // Bookkeep MSSM neutralinos as leptons
    if (abs(hardOutgoing2[i]) == 1000022) nFin++;
    // Bookkeep sleptons as leptons
    if ( abs(hardOutgoing2[i]) == 1000011 || abs(hardOutgoing2[i]) == 2000011
      || abs(hardOutgoing2[i]) == 1000013 || abs(hardOutgoing2[i]) == 2000013
      || abs(hardOutgoing2[i]) == 1000015 || abs(hardOutgoing2[i]) == 2000015)
      nFin++;
  }

  // For very loose hard process definition, check number of hard process
  // lepton explicitly.
  // Check lepton / neutrino containers as leptons
  for(int i =0; i< int(hardOutgoing1.size()); ++i)
    if (hardOutgoing1[i] == 1100)
      for(int j =0; j< int(PosOutgoing1.size()); ++j)
        if (  abs(state[PosOutgoing1[j]].id()) == 11
          || abs(state[PosOutgoing1[j]].id()) == 13
          || abs(state[PosOutgoing1[j]].id()) == 15 )
          nFin++;
  for(int i =0; i< int(hardOutgoing2.size()); ++i)
    if (hardOutgoing2[i] == 1200)
      for(int j =0; j< int(PosOutgoing2.size()); ++j)
        if (  abs(state[PosOutgoing2[j]].id()) == 12
          || abs(state[PosOutgoing2[j]].id()) == 14
          || abs(state[PosOutgoing2[j]].id()) == 16 )
          nFin++;
  return nFin;
}

//--------------------------------------------------------------------------

// Function to get the number of uncoloured final state particles in the
// hard process

int HardProcess::nBosonsOut(){
  int nFin =0;
  for(int i =0; i< int(hardOutgoing1.size()); ++i){
    if (abs(hardOutgoing1[i]) > 20 && abs(hardOutgoing1[i]) <= 25) nFin++;
  }
  for(int i =0; i< int(hardOutgoing2.size()); ++i){
    if (abs(hardOutgoing2[i]) > 20 && abs(hardOutgoing2[i]) <= 25) nFin++;
    if ( hardOutgoing2[i] == 2400) nFin++;
  }
  return nFin;
}

//--------------------------------------------------------------------------

// Function to get the number of coloured initial state partons in the
// hard process

int HardProcess::nQuarksIn(){
  int nIn =0;
  if (hardIncoming1 == 2212 || abs(hardIncoming1) < 10) nIn++;
  if (hardIncoming2 == 2212 || abs(hardIncoming2) < 10) nIn++;
  return nIn;
}

//--------------------------------------------------------------------------

// Function to get the number of uncoloured initial state particles in the
// hard process

int HardProcess::nLeptonIn(){
  int nIn =0;
  if (abs(hardIncoming1) > 10 && abs(hardIncoming1) < 20) nIn++;
  if (abs(hardIncoming2) > 10 && abs(hardIncoming2) < 20) nIn++;
  return nIn;
}

//--------------------------------------------------------------------------

// Function to report if a resonace decay was found in the
// 2->2 hard sub-process in the current state

int HardProcess::hasResInCurrent(){
  for(int i =0; i< int(PosIntermediate.size()); ++i)
    if (PosIntermediate[i] == 0) return 0;
  // Do not count final state bosons as resonaces
  for(int i =0; i< int(PosIntermediate.size()); ++i){
    for(int j =0; j< int(PosOutgoing1.size()); ++j){
      if ( PosIntermediate[i] == PosOutgoing1[j] )
        return 0;
    }
    for(int j =0; j< int(PosOutgoing2.size()); ++j){
      if ( PosIntermediate[i] == PosOutgoing2[j] )
        return 0;
    }
  }
  return 1;
}

//--------------------------------------------------------------------------

// Function to report the number of resonace decays in the 2->2 sub-process
// of the  current state

int HardProcess::nResInCurrent(){
  int nRes = 0;
  for(int i =0; i< int(PosIntermediate.size()); ++i){
    if (PosIntermediate[i] != 0) {
      bool matchesFinalBoson = false;
      for(int j =0; j< int(PosOutgoing1.size()); ++j){
        if ( PosIntermediate[i] == PosOutgoing1[j] )
          matchesFinalBoson = true;
      }
      for(int j =0; j< int(PosOutgoing2.size()); ++j){
        if ( PosIntermediate[i] == PosOutgoing2[j] )
          matchesFinalBoson = true;
      }
      if (!matchesFinalBoson) nRes++;
    }
  }
  return nRes;
}

//--------------------------------------------------------------------------

// Function to report if a resonace decay was found in the
// 2->2 hard core process

bool HardProcess::hasResInProc(){

  for(int i =0; i< int(hardIntermediate.size()); ++i)
    if (hardIntermediate[i] == 0) return false;
  // Do not count final state bosons as resonaces
  for(int i =0; i< int(hardIntermediate.size()); ++i){
    for(int j =0; j< int(hardOutgoing1.size()); ++j){
      if ( hardIntermediate[i] == hardOutgoing1[j] )
        return false;
    }
    for(int j =0; j< int(hardOutgoing2.size()); ++j){
      if ( hardIntermediate[i] == hardOutgoing2[j] )
        return false;
    }
  }
  return true;
}

//--------------------------------------------------------------------------

// Function to print the hard process (for debug)

void HardProcess::list() const {
  cout << "   Hard Process: ";
  cout << " \t " << hardIncoming1 << " + " << hardIncoming2;
  cout << " \t -----> \t ";
  for(int i =0; i < int(hardIntermediate.size());++i)
    cout << hardIntermediate[i] << " ";
  cout << " \t -----> \t ";
  for(int i =0; i < int(hardOutgoing1.size());++i)
    cout << hardOutgoing1[i] << " ";
  for(int i =0; i < int(hardOutgoing2.size());++i)
    cout << hardOutgoing2[i] << " ";
  cout << endl;
}

//--------------------------------------------------------------------------

// Function to list the hard process candiates in the matrix element state
// (for debug)

void HardProcess::listCandidates() const {
  cout << "   Hard Process candidates: ";
  cout << " \t " << hardIncoming1 << " + " << hardIncoming2;
  cout << " \t -----> \t ";
  for(int i =0; i < int(PosIntermediate.size());++i)
    cout << PosIntermediate[i] << " ";
  cout << " \t -----> \t ";
  for(int i =0; i < int(PosOutgoing1.size());++i)
    cout << PosOutgoing1[i] << " ";
  for(int i =0; i < int(PosOutgoing2.size());++i)
    cout << PosOutgoing2[i] << " ";
  cout << endl;
}

//--------------------------------------------------------------------------

// Function to clear hard process information

void HardProcess::clear() {

  // Clear flavour of the first incoming particle
  hardIncoming1 = hardIncoming2 = 0;
  // Clear outgoing particles
  hardOutgoing1.resize(0);
  hardOutgoing2.resize(0);
  // Clear intermediate bosons in the hard 2->2
  hardIntermediate.resize(0);

  // Clear reference event
  state.clear();

  // Clear potential positions of outgoing particles in reference event
  PosOutgoing1.resize(0);
  PosOutgoing2.resize(0);
  // Clear potential positions of intermediate bosons in reference event
  PosIntermediate.resize(0);

  // Clear merging scale read from LHE file
  tms = 0.;

}

//==========================================================================

// The MergingHooks class.

//--------------------------------------------------------------------------

// Initialise MergingHooks class

void MergingHooks::init( Settings settings, Info* infoPtrIn,
  ParticleData* particleDataPtrIn, PartonSystems* partonSystemsPtrIn,
  ostream& os){

  // Save pointers
  infoPtr               = infoPtrIn;
  particleDataPtr       = particleDataPtrIn;
  partonSystemsPtr      = partonSystemsPtrIn;

  // Initialise AlphaS objects for reweighting
  double alphaSvalueFSR = settings.parm("TimeShower:alphaSvalue");
  int    alphaSorderFSR = settings.mode("TimeShower:alphaSorder");
  int    alphaSnfmax    = settings.mode("StandardModel:alphaSnfmax");
  int    alphaSuseCMWFSR= settings.flag("TimeShower:alphaSuseCMW");
  AlphaS_FSRSave.init(alphaSvalueFSR, alphaSorderFSR, alphaSnfmax,
    alphaSuseCMWFSR);
  double alphaSvalueISR = settings.parm("SpaceShower:alphaSvalue");
  int    alphaSorderISR = settings.mode("SpaceShower:alphaSorder");
  int    alphaSuseCMWISR= settings.flag("SpaceShower:alphaSuseCMW");
  AlphaS_ISRSave.init(alphaSvalueISR, alphaSorderISR, alphaSnfmax,
    alphaSuseCMWISR);

  // Initialise AlphaEM objects for reweighting
  int    alphaEMFSRorder = settings.mode("TimeShower:alphaEMorder");
  AlphaEM_FSRSave.init(alphaEMFSRorder, &settings);

  // Initialise merging switches
  doUserMergingSave      = settings.flag("Merging:doUserMerging");
  // Initialise automated MadGraph kT merging
  doMGMergingSave        = settings.flag("Merging:doMGMerging");
  // Initialise kT merging
  doKTMergingSave        = settings.flag("Merging:doKTMerging");
  // Initialise evolution-pT merging
  doPTLundMergingSave    = settings.flag("Merging:doPTLundMerging");
  // Initialise \Delta_R_{ij}, pT_i Q_{ij} merging
  doCutBasedMergingSave  = settings.flag("Merging:doCutBasedMerging");
  // Initialise exact definition of kT
  ktTypeSave             = settings.mode("Merging:ktType");

  // Initialise NL3 switches.
  doNL3TreeSave          = settings.flag("Merging:doNL3Tree");
  doNL3LoopSave          = settings.flag("Merging:doNL3Loop");
  doNL3SubtSave          = settings.flag("Merging:doNL3Subt");
  bool doNL3             = doNL3TreeSave || doNL3LoopSave || doNL3SubtSave;

  // Initialise UNLOPS switches.
  doUNLOPSTreeSave      =  settings.flag("Merging:doUNLOPSTree");
  doUNLOPSLoopSave      =  settings.flag("Merging:doUNLOPSLoop");
  doUNLOPSSubtSave      =  settings.flag("Merging:doUNLOPSSubt");
  doUNLOPSSubtNLOSave   =  settings.flag("Merging:doUNLOPSSubtNLO");
  bool doUNLOPS         = doUNLOPSTreeSave || doUNLOPSLoopSave
                       || doUNLOPSSubtSave || doUNLOPSSubtNLOSave;

  // Initialise UMEPS switches
  doUMEPSTreeSave      =  settings.flag("Merging:doUMEPSTree");
  doUMEPSSubtSave      =  settings.flag("Merging:doUMEPSSubt");
  nReclusterSave       =  settings.mode("Merging:nRecluster");
  nQuarksMergeSave     =  settings.mode("Merging:nQuarksMerge");
  bool doUMEPS         =  doUMEPSTreeSave || doUMEPSSubtSave;

  // Flag to only do phase space cut.
  doEstimateXSection   =  settings.flag("Merging:doXSectionEstimate");

  // Get core process from user input
  processSave           = settings.word("Merging:Process");

  // Clear hard process
  hardProcess.clear();

  // Initialise input event.
  inputEvent.init("(hard process)", particleDataPtr);
  doRemoveDecayProducts = settings.flag("Merging:mayRemoveDecayProducts");

  // Initialise the hard process
  if ( doMGMergingSave )
    hardProcess.initOnLHEF(lheInputFile, particleDataPtr);
  else
    hardProcess.initOnProcess(processSave, particleDataPtr);

  // Parameters for reconstruction of evolution scales
  includeMassiveSave        = settings.flag("Merging:includeMassive");
  enforceStrongOrderingSave = settings.flag("Merging:enforceStrongOrdering");
  scaleSeparationFactorSave = settings.parm("Merging:scaleSeparationFactor");
  orderInRapiditySave       = settings.flag("Merging:orderInRapidity");

  // Parameters for choosing history probabilistically
  nonJoinedNormSave     = settings.parm("Merging:nonJoinedNorm");
  fsrInRecNormSave      = settings.parm("Merging:fsrInRecNorm");
  pickByFullPSave       = settings.flag("Merging:pickByFullP");
  pickByPoPT2Save       = settings.flag("Merging:pickByPoPT2");
  includeRedundantSave  = settings.flag("Merging:includeRedundant");

  // Parameters for scale choices
  unorderedScalePrescipSave   =
    settings.mode("Merging:unorderedScalePrescrip");
  unorderedASscalePrescipSave =
    settings.mode("Merging:unorderedASscalePrescrip");
  unorderedPDFscalePrescipSave =
    settings.mode("Merging:unorderedPDFscalePrescrip");
  incompleteScalePrescipSave  =
    settings.mode("Merging:incompleteScalePrescrip");

  // Parameter for allowing swapping of one colour index while reclustering
  allowColourShufflingSave  =
    settings.flag("Merging:allowColourShuffling");

  // Parameters to allow setting hard process scales to default (dynamical)
  // Pythia values.
  resetHardQRenSave     =  settings.flag("Merging:usePythiaQRenHard");
  resetHardQFacSave     =  settings.flag("Merging:usePythiaQFacHard");

  // Parameters for choosing history by sum(|pT|)
  pickBySumPTSave       = settings.flag("Merging:pickBySumPT");
  herwigAcollFSRSave    = settings.parm("Merging:aCollFSR");
  herwigAcollISRSave    = settings.parm("Merging:aCollISR");

  // Information on the shower cut-off scale
  pT0ISRSave            = settings.parm("SpaceShower:pT0Ref");
  pTcutSave             = settings.parm("SpaceShower:pTmin");
  pTcutSave             = max(pTcutSave,pT0ISRSave);

  // Initialise CKKWL weight
  weightCKKWLSave = 1.;
  weightFIRSTSave = 0.;
  nMinMPISave = 100;
  muMISave = -1.;

  // Initialise merging scale
  tmsValueSave = 0.;
  tmsListSave.resize(0);

  kFactor0jSave         = settings.parm("Merging:kFactor0j");
  kFactor1jSave         = settings.parm("Merging:kFactor1j");
  kFactor2jSave         = settings.parm("Merging:kFactor2j");

  muFSave               = settings.parm("Merging:muFac");
  muRSave               = settings.parm("Merging:muRen");
  muFinMESave           = settings.parm("Merging:muFacInME");
  muRinMESave           = settings.parm("Merging:muRenInME");

  doWClusteringSave     = settings.flag("Merging:allowWClustering");
  doSQCDClusteringSave  = settings.flag("Merging:allowSQCDClustering");
  DparameterSave        = settings.parm("Merging:Dparameter");

  // Save merging scale on maximal number of jets
  if (  doKTMergingSave || doUserMergingSave || doPTLundMergingSave
    || doUMEPS ) {
    // Read merging scale (defined in kT) from input parameter.
    tmsValueSave    = settings.parm("Merging:TMS");
    nJetMaxSave     = settings.mode("Merging:nJetMax");
    nJetMaxNLOSave  = -1;
  } else if (doMGMergingSave) {
    // Read merging scale (defined in kT) from LHE file.
    tmsValueSave    = hardProcess.tms;
    nJetMaxSave     = settings.mode("Merging:nJetMax");
    nJetMaxNLOSave  = -1;
  } else if (doCutBasedMergingSave) {

    // Save list of cuts defining the merging scale.
    nJetMaxSave     = settings.mode("Merging:nJetMax");
    nJetMaxNLOSave  = -1;
    // Write tms cut values to list of cut values,
    // ordered by DeltaR_{ij}, pT_{i}, Q_{ij}.
    tmsListSave.resize(0);
    double drms     = settings.parm("Merging:dRijMS");
    double ptms     = settings.parm("Merging:pTiMS");
    double qms      = settings.parm("Merging:QijMS");
    tmsListSave.push_back(drms);
    tmsListSave.push_back(ptms);
    tmsListSave.push_back(qms);

  }

  if ( doNL3 || doUNLOPS || doEstimateXSection ) {
    tmsValueSave    = settings.parm("Merging:TMS");
    nJetMaxSave     = settings.mode("Merging:nJetMax");
    nJetMaxNLOSave  = settings.mode("Merging:nJetMaxNLO");
  }

  bool writeBanner =  doKTMergingSave || doMGMergingSave
                   || doUserMergingSave
                   || doNL3 || doUNLOPS || doUMEPS
                   || doPTLundMergingSave || doCutBasedMergingSave;

  if (!writeBanner) return;

  // Write banner
  os << "\n *------------------ MEPS Merging Initialization  ---------------"
     << "---*";
  os << "\n |                                                               "
     << "   |\n";
  if ( doKTMergingSave || doMGMergingSave || doUserMergingSave
    || doPTLundMergingSave || doCutBasedMergingSave )
    os << " | CKKW-L merge                                                  "
       << "   |\n"
       << " |"<< setw(34) << processSave << "  with up to"
       << setw(3) << nJetMaxSave << " additional jets |\n";
  else if ( doNL3 )
    os << " | NL3 merge                                                     "
       << "   |\n"
       << " |" << setw(31) << processSave << " with jets up to"
       << setw(3) << nJetMaxNLOSave << " correct to NLO |\n"
       << " | and up to" << setw(3) << nJetMaxSave
       << " additional jets included by CKKW-L merging at LO    |\n";
  else if ( doUNLOPS )
    os << " | UNLOPS merge                                                  "
       << "   |\n"
       << " |" << setw(31) << processSave << " with jets up to"
       << setw(3)<< nJetMaxNLOSave << " correct to NLO |\n"
       << " | and up to" << setw(3) << nJetMaxSave
       << " additional jets included by UMEPS merging at LO     |\n";
  else if ( doUMEPS )
    os << " | UMEPS merge                                                   "
       << "   |\n"
       << " |" << setw(34) << processSave << "  with up to"
       << setw(3) << nJetMaxSave << " additional jets |\n";

  if ( doKTMergingSave )
    os << " | Merging scale is defined in kT, with value ktMS = "
       << tmsValueSave << " GeV";
  else if ( doMGMergingSave )
    os << " | Perform automanted MG/ME merging \n"
       << " | Merging scale is defined in kT, with value ktMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";
  else if ( doUserMergingSave )
    os << " | Merging scale is defined by the user, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << "     |";
  else if ( doPTLundMergingSave )
    os << " | Merging scale is defined by Lund pT, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";
  else if ( doCutBasedMergingSave )
    os << " | Merging scale is defined by combination of Delta R_{ij}, pT_i "
       << "   |\n"
       << " | and Q_{ij} cut, with values                                   "
       << "   |\n"
       << " | Delta R_{ij,min} = "
       << setw(7) << scientific << setprecision(2) << tmsListSave[0]
       << "                                      |\n"
       << " | pT_{i,min}       = "
       << setw(6) << fixed << setprecision(1) << tmsListSave[1]
       << " GeV                                    |\n"
       << " | Q_{ij,min}       = "
       << setw(6) << fixed << setprecision(1) << tmsListSave[2]
       << " GeV                                    |";
  else if ( doNL3TreeSave )
    os << " | Generate tree-level O(alpha_s)-subtracted events              "
       << "   |\n"
       << " | Merging scale is defined by Lund pT, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";
  else if ( doNL3LoopSave )
    os << " | Generate virtual correction unit-weight events                "
       << "   |\n"
       << " | Merging scale is defined by Lund pT, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";
  else if ( doNL3SubtSave )
    os << " | Generate reclustered tree-level events                        "
       << "   |\n"
       << " | Merging scale is defined by Lund pT, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";
  else if ( doUNLOPSTreeSave )
    os << " | Generate tree-level O(alpha_s)-subtracted events              "
       << "   |\n"
       << " | Merging scale is defined by Lund pT, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";
  else if ( doUNLOPSLoopSave )
    os << " | Generate virtual correction unit-weight events                "
       << "   |\n"
       << " | Merging scale is defined by Lund pT, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";
  else if ( doUNLOPSSubtSave )
    os << " | Generate reclustered tree-level events                        "
       << "   |\n"
       << " | Merging scale is defined by Lund pT, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";
  else if ( doUNLOPSSubtNLOSave )
    os << " | Generate reclustered loop-level events                        "
       << "   |\n"
       << " | Merging scale is defined by Lund pT, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";
  else if ( doUMEPSTreeSave )
    os << " | Generate tree-level events                                    "
       << "   |\n"
       << " | Merging scale is defined by Lund pT, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";
  else if ( doUMEPSSubtSave )
    os << " | Generate reclustered tree-level events                        "
       << "   |\n"
       << " | Merging scale is defined by Lund pT, with value tMS = "
       << setw(6) << fixed << setprecision(1) << tmsValueSave << " GeV |";

  os << "\n |                                                               "
     << "   |";
  os << "\n *-------------- END MEPS Merging Initialization  ---------------"
     << "---*\n\n";

}

//--------------------------------------------------------------------------

// Function to check if emission should be rejected.

bool MergingHooks::doVetoEmission( const Event& event) {

  // Do nothing in trial showers, or after first step.
  if ( doIgnoreEmissionsSave ) return false;

  // Do nothing in CKKW-L
  if (  doUserMerging() || doMGMerging() || doKTMerging()
    ||  doPTLundMerging() || doCutBasedMerging() )
     return false;

  // For NLO merging, count and veto emissions above the merging scale
  bool veto = false;
  // Get number of clustering steps
  int nSteps  = getNumberOfClusteringSteps(event);
  // Get merging scale in current event
  double tnow = rhoms( event, false);

  // Get maximal number of additional jets
  int nJetMax = nMaxJets();
  // Always remove emissions above the merging scale for
  // samples containing reclusterings!
  if ( nRecluster() > 0 ) nSteps = max(1, min(nJetMax-2, 1));
  // Check veto condition
  if ( nSteps - 1 < nJetMax && nSteps >= 1 && tnow > tms() ) veto = true;

  // Do not veto if state already includes MPI.
  if ( infoPtr->nMPI() > 1 ) veto = false;

  // When performing NL3 merging of tree-level events, reset the
  // CKKWL weight.
  if ( veto && doNL3Tree() ) setWeightCKKWL(0.);

  // If the emission is allowed, do not check any further emissions
  if ( !veto ) doIgnoreEmissionsSave = true;

  // Done
  return veto;

}

//--------------------------------------------------------------------------

// Function to check if emission should be rejected.

bool MergingHooks::doVetoStep( const Event& process, const Event& event,
  bool doResonance ) {

  // Do nothing in trial showers, or after first step.
  if ( doIgnoreStepSave && !doResonance ) return false;

  // Do nothing in UMEPS or UNLOPS
  if ( doUMEPSTree() || doUMEPSSubt() || doUMEPSMerging() || doUNLOPSTree()
    || doUNLOPSLoop() || doUNLOPSSubt() || doUNLOPSSubtNLO()
    || doUNLOPSMerging() )
    return false;

  // Get number of clustering steps. If necessary, remove resonance
  // decay products first.
  int nSteps  = (doResonance) ? getNumberOfClusteringSteps(process)
              : getNumberOfClusteringSteps( bareEvent( process, false) );
  // Get maximal number of additional jets.
  int nJetMax = nMaxJets();
  // Get merging scale in current event.
  double tnow = tmsNow( event );

  // For non-resonant showers, simply check veto. If event should indeed be
  // vetoed, save the current pT and weights in case the veto needs to be
  // revoked.
  if ( !doResonance ) {

    // Store pT to check if veto needs to be revoked later.
    pTsave = infoPtr->pTnow();
    if ( nRecluster() == 1) nSteps--;

    // Check merging veto condition.
    bool veto = false;
    if ( nSteps > nMaxJetsNLO() && nSteps < nJetMax && tnow > tms() ) {
      // Set weight to zero if event should be vetoed.
      weightCKKWL1Save = 0.;
      // Save weight before veto, in case veto needs to be revoked.
      weightCKKWL2Save = getWeightCKKWL();
      // Reset stored weights.
      setWeightCKKWL(0.);
      veto = true;
    }

    // Done
    return veto;

  // For resonant showers, check if any previous veto should be revoked.
  // This means we treat showers off resonance decay products identical to
  // MPI: If a hard resonance emission has been produced, the event should
  // have been kept. Following this reasoning, it might be necessary to revoke
  // any previous veto.
  } else {

    // Initialise switch to revoke vetoing.
    bool revokeVeto = false;

    // Nothing to check if pTsave was not stored, i.e. no emission to
    // possibly veto was recorded.
    // Only allow revoking the veto for diboson processes, with resonant
    // electroweak bosons
    bool check =  (nHardInLeptons() == 0)&& (nHardOutLeptons() == 2)
               && (nHardOutPartons() == 2);

    // For current purpose only!!!
    check = false;

    // For hadronic resonance decays at hadron colliders, do not veto
    // events with a hard emission of the resonance decay products,
    // since such states are not included in the matrix element
    if ( pTsave > 0. && check ) {

      // Check how many resonance decay systems are allowed
      int nResNow = nResInCurrent();

      // Find systems really containing only emission off a resonance
      // decay
      vector<int>goodSys;
      // Resonance decay systems are considered last, thus at the end of
      // the list
      int sysSize = partonSystemsPtr->sizeSys();
      for ( int i=0; i < nResNow; ++i ) {
        if ( partonSystemsPtr->sizeOut(sysSize - 1 - i) == 3 )
          goodSys.push_back(sysSize - 1 - i);
      }

      // Check the members of the resonance decay systems
      for ( int i=0; i < int(goodSys.size()); ++i ) {

        // Save the (three) members of the resonance decay system
        int iMem1 = partonSystemsPtr->getOut(goodSys[i],0);
        int iMem2 = partonSystemsPtr->getOut(goodSys[i],1);
        int iMem3 = partonSystemsPtr->getOut(goodSys[i],2);

        // Now find emitted gluon or gamma
        int iEmtGlue = ((event[iMem1].id() == 21) ? iMem1
                     : ((event[iMem2].id() == 21) ? iMem2
                       : ((event[iMem3].id() == 21) ? iMem3: 0)));
        int iEmtGamm = ((event[iMem1].id() == 22) ? iMem1
                     : ((event[iMem2].id() == 22) ? iMem2
                       : ((event[iMem3].id() == 22) ? iMem3: 0)));
        // Per system, only one emission
        int iEmt = (iEmtGlue != 0) ? iEmtGlue : iEmtGamm;

        int iRad = 0;
        int iRec = 0;
        if ( iEmt == iMem1 ) {
          iRad = (event[iMem2].mother1() != event[iMem2].mother2())
               ? iMem2 : iMem3;
          iRec = (event[iMem3].mother1() == event[iMem3].mother2())
               ? iMem3 : iMem2;
        } else if ( iEmt == iMem2 ) {
          iRad = (event[iMem1].mother1() != event[iMem1].mother2())
               ? iMem1 : iMem3;
          iRec = (event[iMem3].mother1() == event[iMem3].mother2())
               ? iMem3 : iMem1;
        } else {
          iRad = (event[iMem1].mother1() != event[iMem1].mother2())
               ? iMem1 : iMem2;
          iRec = (event[iMem2].mother1() == event[iMem2].mother2())
               ? iMem2 : iMem1;
        }

        double pTres = rhoPythia(event[iRad], event[iEmt], event[iRec], 1);

        // Revoke previous veto of last emission if a splitting of the
        // resonance produced a harder parton, i.e. we are inside the
        // PS region
        if ( pTres > pTsave ) {
          revokeVeto = true;
        // Do nothing (i.e. allow other first emission veto) for soft
        // splitting
        } else {
          revokeVeto = false;
        }
      // Done with one system
      }
    // Done with all systems
    }

    // Check veto condition
    bool veto = false;
    if ( revokeVeto && check ) {
      setWeightCKKWL(weightCKKWL2Save);
    } else if ( check ) {
      setWeightCKKWL(weightCKKWL1Save);
      if ( weightCKKWL1Save == 0. ) veto = true;
    }

    // Check veto condition.
    if ( !check && nSteps > nMaxJetsNLO() && nSteps < nJetMax && tnow > tms()){
      // Set stored weights to zero.
      setWeightCKKWL(0.);
      // Now allow veto.
      veto = true;
    }

    // If the emission is allowed, do not check any further emissions
    if ( !veto || !doIgnoreStepSave ) doIgnoreStepSave = true;

    // Done
    return veto;

  }

  // Done
  return false;

}

//--------------------------------------------------------------------------

// Return event stripped from decay products.

Event MergingHooks::bareEvent(const Event& inputEventIn,
  bool storeInputEvent ) {

  // Find and detach decay products.
  Event newProcess = Event();
  newProcess.init("(hard process-modified)", particleDataPtr);

  // If desired, store input event.
  if ( storeInputEvent ) {
    resonances.resize(0);
    inputEvent.clear();
    for (int i = 0; i < inputEventIn.size(); ++i)
      inputEvent.append( inputEventIn[i] );
    for (int i = 0; i < inputEventIn.sizeJunction(); ++i)
      inputEvent.appendJunction( inputEventIn.getJunction(i) );
    inputEvent.saveSize();
    inputEvent.saveJunctionSize();
  }

  // Now remove decay products.
  if ( doRemoveDecayProducts ) {

    // Add the beam and initial partons to the event record.
    for (int i = 0; i < inputEventIn.size(); ++ i) {
      if ( inputEventIn[i].mother1() > 4
        || inputEventIn[i].statusAbs() == 22
        || inputEventIn[i].statusAbs() == 23)
        break;
      newProcess.append(inputEventIn[i]);
    }

    // Add the intermediate particles to the event record.
    for (int i = 0; i < inputEventIn.size(); ++ i) {
      if (inputEventIn[i].mother1() > 4) break;
      if ( inputEventIn[i].status() == -22) {
        int j = newProcess.append(inputEventIn[i]);
        newProcess[j].statusPos();
        if ( storeInputEvent ) resonances.push_back( make_pair(j, i) );
        newProcess[j].daughters(0, 0);
      }
    }

    // Add remaining outgoing particles to the event record.
    for (int i = 0; i < inputEventIn.size(); ++ i) {
      if (inputEventIn[i].mother1() > 4) break;
      if ( inputEventIn[i].statusAbs() != 11
        && inputEventIn[i].statusAbs() != 12
        && inputEventIn[i].statusAbs() != 21
        && inputEventIn[i].statusAbs() != 22)
        newProcess.append(inputEventIn[i]);
    }

    // Update event colour tag to maximum in whole process.
    int maxColTag = 0;
    for (int i = 0; i < inputEventIn.size(); ++ i) {
      if ( inputEventIn[i].col() > maxColTag )
        maxColTag = inputEventIn[i].col();
      if ( inputEventIn[i].acol() > maxColTag )
        maxColTag = inputEventIn[i].acol();
    }
    newProcess.initColTag(maxColTag);

    // Copy junctions from process to newProcess.
    for (int i = 0; i < inputEventIn.sizeJunction(); ++i)
      newProcess.appendJunction( inputEventIn.getJunction(i));

    newProcess.saveSize();
    newProcess.saveJunctionSize();

  } else {
    newProcess = inputEventIn;
  }

  // Remember scale
  newProcess.scale( inputEventIn.scale() );

  // Done
  return newProcess;

}

//--------------------------------------------------------------------------

// Write event with decay products attached to argument. Only possible if an
// input event with decay producs had been stored before.

bool MergingHooks::reattachResonanceDecays(Event& process ) {

  // Now reattach the decay products.
  if ( doRemoveDecayProducts && inputEvent.size() > 0 ) {

    int sizeBef = process.size();
    // Vector of resonances for which the decay products were already attached.
    vector<int> iAftChecked;
    // Reset daughters and status of intermediate particles.
    for ( int i = 0; i < int(resonances.size()); ++i ) {
      for (int j = 0; j < sizeBef; ++j ) {
        if ( j != resonances[i].first ) continue;

        int iOldDaughter1 = inputEvent[resonances[i].second].daughter1();
        int iOldDaughter2 = inputEvent[resonances[i].second].daughter2();

        // Get momenta in case of reclustering.
        int iHardMother      = resonances[i].second;
        Particle& hardMother = inputEvent[iHardMother];
        // Find current mother copy (after clustering).
        int iAftMother       = 0;
        for ( int k = 0; k < process.size(); ++k )
          if ( process[k].id() == inputEvent[resonances[i].second].id() ) {
            // Only attempt if decays of this resonance were not attached.
            bool checked = false;
            for ( int l = 0; l < int(iAftChecked.size()); ++l )
              if ( k == iAftChecked[l] )
                checked = true;
            if ( !checked ) {
              iAftChecked.push_back(k);
              iAftMother = k;
              break;
            }
          }

        Particle& aftMother  = process[iAftMother];

        // Resonance can have been moved by clustering,
        // so prepare to update colour and momentum information for system.
        int colBef  = hardMother.col();
        int acolBef = hardMother.acol();
        int colAft  = aftMother.col();
        int acolAft = aftMother.acol();
        RotBstMatrix M;
        M.bst( hardMother.p(), aftMother.p());

        // Attach resonance decay products.
        int iNewDaughter1 = 0;
        int iNewDaughter2 = 0;
        for ( int k = iOldDaughter1; k <= iOldDaughter2; ++k ) {
          if ( k == iOldDaughter1 )
            iNewDaughter1 = process.append(inputEvent[k] );
          else
            iNewDaughter2 = process.append(inputEvent[k] );
          process.back().statusPos();
          Particle& now = process.back();
          // Update colour and momentum information.
          if (now.col()  != 0 && now.col()  == colBef ) now.col(colAft);
          if (now.acol() != 0 && now.acol() == acolBef) now.acol(acolAft);
          now.rotbst( M);
          // Update vertex information.
          if (now.hasVertex()) now.vProd( aftMother.vDec() );
          // Update mothers.
          now.mothers(iAftMother,0);
        }

        process[iAftMother].daughters( iNewDaughter1, iNewDaughter2 );
        process[iAftMother].statusNeg();

        // Loop through event and attach remaining decays
        int iDec = 0;
        do {
          if ( process[iDec].isFinal() && process[iDec].canDecay()
            && process[iDec].mayDecay() && process[iDec].isResonance() ) {

            int iD1 = process[iDec].daughter1();
            int iD2 = process[iDec].daughter2();

            // Done if no daughters exist.
            if ( iD1 == 0 || iD2 == 0 ) continue;

            // Attach daughters.
            int iNewDaughter12 = 0;
            int iNewDaughter22 = 0;
            for ( int k = iD1; k <= iD2; ++k ) {
              if ( k == iD1 )
                iNewDaughter12 = process.append(inputEvent[k] );
              else
                iNewDaughter22 = process.append(inputEvent[k] );
              process.back().statusPos();
              Particle& now = process.back();
              // Update colour and momentum information.
              if (now.col() != 0 && now.col() == colBef ) now.col(colAft);
              if (now.acol()!= 0 && now.acol()== acolBef) now.acol(acolAft);
              now.rotbst( M);
              // Update vertex information.
              if (now.hasVertex()) now.vProd( process[iDec].vDec() );
              // Update mothers.
              now.mothers(iDec,0);
            }

            // Modify mother status and daughters.
            process[iDec].status(-22);
            process[iDec].daughters(iNewDaughter12, iNewDaughter22);

          // End of loop over all entries.
          }
        } while (++iDec < process.size());
      } // End loop over process entries.
    } // End loop over resonances.

    // Update event colour tag to maximum in whole process.
    int maxColTag = 0;
    for (int i = 0; i < process.size(); ++ i) {
      if (process[i].col() > maxColTag) maxColTag = process[i].col();
      if (process[i].acol() > maxColTag) maxColTag = process[i].acol();
    }
    process.initColTag(maxColTag);

  }

  // Done.
  return (doRemoveDecayProducts) ? inputEvent.size() > 0 : true;

}

//--------------------------------------------------------------------------

bool MergingHooks::isInHard( int iPos, const Event& event){

  // MPI not part of hard process
  if ( event[iPos].statusAbs() > 30 && event[iPos].statusAbs() < 40 )
    return false;
  // Beam remnants and hadronisation not part of hard process
  if ( event[iPos].statusAbs() > 60 )
    return false;

  // Still MPI: Check that the particle is not due to radiation off MPIs.
  // First get all intermediate MPI partons in the state.
  vector<int> mpiParticlePos;
  mpiParticlePos.clear();
  for ( int i=0; i < event.size(); ++i )
    if ( event[i].statusAbs() > 30
      && event[i].statusAbs() < 40 )
      mpiParticlePos.push_back(i);
  // Disregard any parton iPos that has MPI ancestors.
  for ( int i=0; i < int(mpiParticlePos.size()); ++i)
    if ( event[iPos].isAncestor( mpiParticlePos[i]) )
      return false;

  // Disallow other systems.
  // Get sub-system of particles for iPos
  int iSys = partonSystemsPtr->getSystemOf(iPos, !event[iPos].isFinal() );
  if ( iSys > 0 ) {
    // Check all partons belonging to the same system as iPos. If any is
    // produced in MPI or has MPI ancestors, the whole system is not the
    // hard subprocess, i.e. iPos is not in the hard subprocess.
    int sysSize = partonSystemsPtr->sizeAll(iSys);
    for ( int i = 0; i < sysSize; ++i ) {
      int iPosNow = partonSystemsPtr->getAll( iSys, i );
      // MPI not part of hard process
      if ( event[iPosNow].statusAbs() > 30
        && event[iPosNow].statusAbs() < 40 )
        return false;
      // Disregard any parton iPos that has MPI ancestors.
      for ( int j=0; j < int(mpiParticlePos.size()); ++j)
        if ( event[iPosNow].isAncestor( mpiParticlePos[j]) )
          return false;
      // Beam remnants and hadronisation not part of hard process
      if ( event[iPosNow].statusAbs() > 60 )
        return false;
    }
  }

  // Check if any ancestor contains the hard incoming partons as daughters.
  // Begin loop to trace upwards from the daughter.
  bool containsInitialParton = false;
  int iUp = iPos;
  for ( ; ; ) {
    // If out of range then failed to find match.
    if (iUp <= 0 || iUp > event.size()) break;
    // If positive match then done.
    if ( iUp == 3 || iUp == 4 ) {
      containsInitialParton = true;
      break;
    }
    if ( event[iUp].mother1() == 1
      && (event[iUp].daughter1() == 3 || event[iUp].daughter2() == 3) ) {
      containsInitialParton = true;
      break;
    }
    if ( event[iUp].mother1() == 2
      && (event[iUp].daughter1() == 4 || event[iUp].daughter2() == 4) ) {
      containsInitialParton = true;
      break;
    }
    // If unique mother then keep on moving up the chain.
    iUp = event[iUp].mother1();
  }

  if ( !containsInitialParton ) return false;

  // Done
  return true;

}

//--------------------------------------------------------------------------

// Function to return the number of clustering steps for the current event

int MergingHooks::getNumberOfClusteringSteps(const Event& event){

  // Count the number of final state partons
  int nFinalPartons = 0;
  for ( int i=0; i < event.size(); ++i)
    if ( event[i].isFinal()
      && isInHard( i, event)
      && (event[i].isQuark() || event[i].isGluon()) )
      nFinalPartons++;

  // Count the number of final state leptons
  int nFinalLeptons = 0;
  for( int i=0; i < event.size(); ++i)
    if ( event[i].isFinal() && isInHard( i, event) && event[i].isLepton())
      nFinalLeptons++;

  // Add neutralinos to number of leptons
  for( int i=0; i < event.size(); ++i)
    if ( event[i].isFinal() && isInHard( i, event)
       && event[i].idAbs() == 1000022)
      nFinalLeptons++;

  // Add sleptons to number of leptons
  for( int i=0; i < event.size(); ++i)
    if ( event[i].isFinal() && isInHard( i, event)
       && (event[i].idAbs() == 1000011
        || event[i].idAbs() == 2000011
        || event[i].idAbs() == 1000013
        || event[i].idAbs() == 2000013
        || event[i].idAbs() == 1000015
        || event[i].idAbs() == 2000015) )
      nFinalLeptons++;

  // Count the number of final state electroweak bosons
  int nFinalBosons = 0;
  for( int i=0; i < event.size(); ++i)
    if ( event[i].isFinal() && isInHard( i, event)
      && ( event[i].idAbs() == 22
        || event[i].idAbs() == 23
        || event[i].idAbs() == 24
        || event[i].idAbs() == 25 ) )
      nFinalBosons++;

  // Save sum of all final state particles
  int nFinal = nFinalPartons + nFinalLeptons
             + 2*(nFinalBosons - nHardOutBosons() );

  // Return the difference to the core process outgoing particles
  return (nFinal - nHardOutPartons() - nHardOutLeptons() );

}

//--------------------------------------------------------------------------

// Function to check if event contains an emission not present in the hard
// process.

bool MergingHooks::isFirstEmission(const Event& event ) {

  // If the beam remnant treatment or hadronisation has already started, do
  // no veto.
  for ( int i=0; i < event.size(); ++i)
    if ( event[i].statusAbs() > 60 ) return false;

  // Count particle types
  int nFinalQuarks   = 0;
  int nFinalGluons   = 0;
  int nFinalLeptons  = 0;
  int nFinalBosons   = 0;
  int nFinalPhotons  = 0;
  int nFinal         = 0;
  for( int i=0; i < event.size(); ++i) {
    if (event[i].isFinal() && isInHard(i, event) ){
      if ( event[i].spinType() == 2 && event[i].colType() == 0)
        nFinalLeptons++;
      if ( event[i].id()    == 23
        || event[i].idAbs() == 24
        || event[i].id()    == 25)
        nFinalBosons++;
      if (  event[i].id()   == 22)
        nFinalPhotons++;
      if ( event[i].isQuark())
        nFinalQuarks++;
      if ( event[i].isGluon())
        nFinalGluons++;
      if ( !event[i].isDiquark() )
        nFinal++;
    }
  }

  // Return highest value if the event does not contain any final state
  // coloured particles.
  if (nFinalQuarks + nFinalGluons == 0) return false;

  // Use MergingHooks functions to get information on the hard process.
  int nLeptons      = nHardOutLeptons();

  // The state is already in the PS region if the number of leptons had been
  // increased bt QED splittings.
  if (nFinalLeptons > nLeptons) return false;

  // If the mumber of photons if larger than in the hard process, QED
  // radiation has pushed the state into the PS region.
  int nPhotons = 0;
  for(int i =0; i< int(hardProcess.hardOutgoing1.size()); ++i)
    if (hardProcess.hardOutgoing1[i] == 22)
      nPhotons++;
  for(int i =0; i< int(hardProcess.hardOutgoing2.size()); ++i)
    if (hardProcess.hardOutgoing2[i] == 22)
      nPhotons++;
  if (nFinalPhotons > nPhotons) return false;

  // Done
  return true;
}

//--------------------------------------------------------------------------

// Function to set the correct starting scales of the shower

// Set starting scales
bool MergingHooks::setShowerStartingScales( bool isTrial,
  bool doMergeFirstEmm, double& pTscaleIn, const Event& event,
  double& pTmaxFSRIn, bool& limitPTmaxFSRIn,
  double& pTmaxISRIn, bool& limitPTmaxISRIn,
  double& pTmaxMPIIn, bool& limitPTmaxMPIIn ) {

  // MPI treated differently in case of qcd djiet merging, since hard MPI
  // can be misidentified as hard process.
  bool isPureQCD = ( getProcessString().compare("pp>jj") == 0 );
  int nSteps     = getNumberOfClusteringSteps( bareEvent(event, false) );
  double scale   = event.scale();

  bool doRecluster = doUMEPSSubt() || doNL3Subt() || doUNLOPSSubt()
                  || doUNLOPSSubtNLO();

  // Set correct starting scales for trial showers.
  if ( isTrial ) {
    // Reset shower scales.
    pTmaxISRIn = pTmaxFSRIn = scale;
    if ( limitPTmaxISRIn ) pTmaxISRIn = min(scale,muF());
    if ( limitPTmaxFSRIn ) pTmaxFSRIn = min(scale,muF());
    // Reset MPI scale.
    if ( !isPureQCD ) pTmaxMPIIn = scale;
    else  pTmaxMPIIn = infoPtr->eCM();
    // Reset phase space limitation flags
    if ( pTscaleIn < infoPtr->eCM() ) {
      limitPTmaxISRIn = limitPTmaxFSRIn = true;
      if ( !isPureQCD ) limitPTmaxMPIIn = true;
      else limitPTmaxMPIIn = false;
    }
  }

  // Reset starting scales.
  if ( isPureQCD && doMergeFirstEmm ) {
    // Set correct shower starting scales.
    if ( nSteps == 0 && !doRecluster && pTscaleIn < infoPtr->eCM() ) {
      pTmaxMPIIn = infoPtr->eCM();
      limitPTmaxMPIIn = false;
    }
  }

  // Reset starting scales.
  if ( !isPureQCD && doMergeFirstEmm ) {
    // Set correct shower starting scales.
    if ( nSteps == 0 && !doRecluster ) {
      pTscaleIn = infoPtr->eCM();
      limitPTmaxMPIIn = false;
    }
    if ( limitPTmaxISRIn ) pTscaleIn = min(scale,muF());
    if ( limitPTmaxFSRIn ) pTscaleIn = min(scale,muF());
    pTmaxISRIn = pTmaxFSRIn = pTscaleIn;
    // Reset MPI starting scale. Standard treatment
    if ( pTscaleIn < infoPtr->eCM()
      &&  !(nSteps == 0 && !doRecluster
        && (limitPTmaxISRIn || limitPTmaxFSRIn)) ) {
      pTmaxMPIIn = pTscaleIn;
      limitPTmaxMPIIn = true;
    }
    if (doRecluster) {
      pTmaxMPIIn = muMI();
      limitPTmaxMPIIn = true;
    }
  }

  // Done
  return true;

}

//--------------------------------------------------------------------------

// Function to return the value of the merging scale function in the current
// event.

double MergingHooks::tmsNow( const Event& event ) {

  // Get merging scale in current event.
  double tnow = 0.;
  // Use KT/Durham merging scale definition.
  if ( doKTMerging()  || doMGMerging() )
    tnow = kTms(event);
  // Use Lund PT merging scale definition.
  else if ( doPTLundMerging() )
    tnow = rhoms(event, false);
  // Use DeltaR_{ij}, pT_i, Q_{ij} combination merging scale definition.
  else if ( doCutBasedMerging() )
    tnow = cutbasedms(event);
  // Use NLO merging (Lund PT) merging scale definition.
  else if ( doNL3Merging() )
    tnow = rhoms(event, false);
  // Use NLO merging (Lund PT) merging scale definition.
  else if ( doUNLOPSMerging() )
    tnow = rhoms(event, false);
  // Use UMEPS (Lund PT) merging scale definition.
  else if ( doUMEPSMerging() )
    tnow = rhoms(event, false);
  // Use user-defined merging scale.
  else
    tnow = tmsDefinition(event);
  // Return merging scale value. Done
  return tnow;
}

//--------------------------------------------------------------------------

// Function to check if the properties of the input particle should be
// checked against the cut-based merging scale defintion.

bool MergingHooks::checkAgainstCut( const Particle& particle){

  // Do not check uncoloured particles.
  if (particle.colType() == 0) return false;
  // By default, use u-, d-, c-, s- and b-quarks.
  if ( particle.idAbs() != 21 && particle.idAbs() > nQuarksMergeSave )
    return false;
  // Done
  return true;

}

//--------------------------------------------------------------------------

// Function to return the minimal kT in the event. If doKTMerging = true, this
// function will be used as a merging scale definition.

double MergingHooks::kTms(const Event& event) {

  // Only check first emission.
  if (!isFirstEmission(event)) return 0.;

  // Find all electroweak decayed bosons in the state.
  vector<int> ewResonancePos;
  ewResonancePos.clear();
  for (int i=0; i < event.size(); ++i)
    if ( abs(event[i].status()) == 22
      && ( event[i].idAbs() == 22
        || event[i].idAbs() == 23
        || event[i].idAbs() == 24
        || event[i].idAbs() == 25
        || event[i].idAbs() == 6 ) )
      ewResonancePos.push_back(i);

  // Declare final parton vectors
  vector <int> FinalPartPos;
  FinalPartPos.clear();
  // Search inEvent record for final state partons.
  // Exclude decay products of ew resonance.
  for (int i=0; i < event.size(); ++i){
    if ( event[i].isFinal()
      && isInHard( i, event )
      && event[i].colType() != 0
      && checkAgainstCut(event[i]) ){
      bool isDecayProduct = false;
      for(int j=0; j < int(ewResonancePos.size()); ++j)
        if ( event[i].isAncestor( ewResonancePos[j]) )
          isDecayProduct = true;
      // Except for e+e- -> jets, do not check radiation in resonance decays.
      if ( !isDecayProduct
        || getProcessString().compare("e+e->jj") == 0
        || getProcessString().compare("e+e->(z>jj)") == 0 )
        FinalPartPos.push_back(i);
    }
  }

  // Find minimal Durham kT in event, using own function: Check
  // definition of separation
  int type = (event[3].colType() == 0
           && event[4].colType() == 0) ? -1 : ktTypeSave;
  // Find minimal kT
  double ktmin = event[0].e();
  for(int i=0; i < int(FinalPartPos.size()); ++i){
    double kt12  = ktmin;
    // Compute separation to the beam axis for hadronic collisions
    if (type == 1 || type == 2) {
      double temp = event[FinalPartPos[i]].pT();
      kt12 = min(kt12, temp);
    }
    // Compute separation to other final state jets
    for(int j=i+1; j < int(FinalPartPos.size()); ++j) {
      double temp = kTdurham( event[FinalPartPos[i]], event[FinalPartPos[j]],
                      type, DparameterSave);
      kt12 = min(kt12, temp);
    }
    // Keep the minimal Durham separation
    ktmin = min(ktmin,kt12);
  }

  // Return minimal Durham kT
  return ktmin;

}

//--------------------------------------------------------------------------

// Function to compute durham y separation from Particle input.

double MergingHooks::kTdurham(const Particle& RadAfterBranch,
  const Particle& EmtAfterBranch, int Type, double D ){

  // Declare return variable
  double ktdur;
  // Save 4-momenta of final state particles
  Vec4 jet1 = RadAfterBranch.p();
  Vec4 jet2 = EmtAfterBranch.p();

  if ( Type == -1) {
    // Get angle between jets for e+e- collisions, make sure that
    // -1 <= cos(theta) <= 1
    double costh;
    if (jet1.pAbs()*jet2.pAbs() <=0.) costh = 1.;
    else {
      costh = costheta(jet1,jet2);
    }
    // Calculate kt durham separation between jets for e+e- collisions
    ktdur = 2.0*min( pow(jet1.e(),2) , (pow(jet2.e(),2)) )*(1.0 - costh);
  } else if ( Type == 1 ){
    // Get delta_y for hadronic collisions:
    // Get mT of first jet
    double mT1sq = jet1.m2Calc() + jet1.pT2();
    double mT1 = 0.;
    if (mT1sq < 0) mT1 = - sqrt(-mT1sq);
    else mT1 = sqrt(mT1sq);
    // Get mT of second jet
    double mT2sq = jet2.m2Calc() + jet2.pT2();
    double mT2 = 0.;
    if (mT2sq < 0) mT2 = - sqrt(-mT2sq);
    else mT2 = sqrt(mT2sq);
    // Get rapidity of first jet
    double y1 = log( ( jet1.e() + abs(jet1.pz()) ) / mT1 );
    if (jet1.pz() < 0) y1 *= -1.;
    // Get rapidity of second jet
    double y2 = log( ( jet2.e() + abs(jet2.pz()) ) / mT2 );
    if (jet2.pz() < 0) y2 *= -1.;
    // Get delta_phi for hadronic collisions
    double pt1 = sqrt( pow(jet1.px(),2) + pow(jet1.py(),2) );
    double pt2 = sqrt( pow(jet2.px(),2) + pow(jet2.py(),2) );
    double cosdPhi = ( jet1.px()*jet2.px() + jet1.py()*jet2.py() ) / (pt1*pt2);
    double dPhi = acos( cosdPhi );
    // Calculate kT durham like fastjet,
    // but with rapidity instead of pseudo-rapidity
    ktdur = min( pow(pt1,2),pow(pt2,2) )
          * ( pow(y1-y2,2) + pow(dPhi,2) ) / pow(D,2);
  } else if ( Type == 2 ){

    // Get mT of first jet
    double mT1sq = jet1.m2Calc() + jet1.pT2();
    double mT1 = 0.;
    if (mT1sq < 0) mT1 = - sqrt(-mT1sq);
    else mT1 = sqrt(mT1sq);
    // Get mT of second jet
    double mT2sq = jet2.m2Calc() + jet2.pT2();
    double mT2 = 0.;
    if (mT2sq < 0) mT2 = - sqrt(-mT2sq);
    else mT2 = sqrt(mT2sq);
    // Get pseudo-rapidity of first jet
    double eta1 = log( ( sqrt(jet1.px()*jet1.px() + jet1.py()*jet1.py()
                            + jet1.pz()*jet1.pz())  + abs(jet1.pz()) ) / mT1);
    if (jet1.pz() < 0) eta1 *= -1.;
    // Get pseudo-rapidity of second jet
    double eta2 = log( ( sqrt(jet2.px()*jet2.px() + jet2.py()*jet2.py()
                            + jet2.pz()*jet2.pz())  + abs(jet2.pz()) ) / mT2);
    if (jet2.pz() < 0) eta2 *= -1.;

    // Get delta_phi and cos(Delta_phi) for hadronic collisions
    double pt1 = sqrt( pow(jet1.px(),2) + pow(jet1.py(),2) );
    double pt2 = sqrt( pow(jet2.px(),2) + pow(jet2.py(),2) );
    double cosdPhi = ( jet1.px()*jet2.px() + jet1.py()*jet2.py() ) / (pt1*pt2);
    double dPhi = acos( cosdPhi );
    // Calculate kT durham like fastjet
    ktdur = min( pow(pt1,2),pow(pt2,2) )
          * ( pow(eta1-eta2,2) + pow(dPhi,2) ) / pow(D,2);
  } else if ( Type == 3 ){
    // Get cosh(Delta_eta) for hadronic collisions
    double eta1 = 0.5*log( (jet1.e() + jet1.pz()) / (jet1.e() - jet1.pz()) );
    double eta2 = 0.5*log( (jet2.e() + jet2.pz()) / (jet2.e() - jet2.pz()) );
    double coshdEta = cosh( eta1 - eta2 );
    // Get delta_phi and cos(Delta_phi) for hadronic collisions
    double pt1 = sqrt( pow(jet1.px(),2) + pow(jet1.py(),2) );
    double pt2 = sqrt( pow(jet2.px(),2) + pow(jet2.py(),2) );
    double cosdPhi = ( jet1.px()*jet2.px() + jet1.py()*jet2.py() ) / (pt1*pt2);
    // Calculate kT durham separation "SHERPA-like"
    ktdur = 2.0*min( pow(pt1,2),pow(pt2,2) )
          * ( coshdEta - cosdPhi ) / pow(D,2);
  } else {
    ktdur = 0.0;
  }
  // Return kT
  return sqrt(ktdur);
}

//--------------------------------------------------------------------------

// Find the minimal Lund pT between coloured partons in the input
// event. If doPTLundMerging = true, this function will be used as a merging
// scale definition.

double MergingHooks::rhoms( const Event& event, bool withColour){

  // Only check first emission.
  if (!isFirstEmission(event)) return 0.;

  // Find all electroweak decayed bosons in the state.
  vector<int> ewResonancePos;
  ewResonancePos.clear();
  for (int i=0; i < event.size(); ++i)
    if ( abs(event[i].status()) == 22
      && ( event[i].idAbs() == 22
        || event[i].idAbs() == 23
        || event[i].idAbs() == 24
        || event[i].idAbs() == 25
        || event[i].idAbs() == 6 ) )
      ewResonancePos.push_back(i);

  // Declare final parton vectors
  vector <int> FinalPartPos;
  FinalPartPos.clear();
  // Search inEvent record for final state partons.
  // Exclude decay products of ew resonance.
  for (int i=0; i < event.size(); ++i){

    if ( event[i].isFinal()
      && isInHard( i, event )
      && event[i].colType() != 0
      && checkAgainstCut(event[i]) ){
      bool isDecayProduct = false;
      for(int j=0; j < int(ewResonancePos.size()); ++j)
        if ( event[i].isAncestor( ewResonancePos[j]) )
          isDecayProduct = true;
      // Except for e+e- -> jets, do not check radiation in resonance decays.
      if ( !isDecayProduct
        || getProcessString().compare("e+e->jj") == 0
        || getProcessString().compare("e+e->(z>jj)") == 0 )
        FinalPartPos.push_back(i);
    }
  }

  // Get index of first incoming
  int in1 = 0;
  for (int i=0; i < event.size(); ++i)
    if (abs(event[i].status()) == 41 ){
      in1 = i;
      break;
    }

  // Get index of second incoming
  int in2 = 0;
  for (int i=0; i < event.size(); ++i)
    if (abs(event[i].status()) == 42 ){
      in2 = i;
      break;
    }

  // If no incoming of the cascade are found, try incoming
  if (in1 == 0 || in2 == 0){
    // Find current incoming partons
    for(int i=3; i < int(event.size()); ++i){
      if ( !isInHard( i, event ) ) continue;
      if (event[i].mother1() == 1) in1 = i;
      if (event[i].mother1() == 2) in2 = i;
    }
  }

  // Find minimal pythia pt in event
  double ptmin = event[0].e();
  for(int i=0; i < int(FinalPartPos.size()); ++i){
    double pt12  = ptmin;
    // Compute pythia ISR separation i-jet and first incoming
    if (event[in1].colType() != 0) {
      double temp = rhoPythia( event[in1],
                      event[FinalPartPos[i]], event[in2], -1 );
      pt12 = min(pt12, temp);
    }

    // Compute pythia ISR separation i-jet and second incoming
    if ( event[in2].colType() != 0) {
      double temp = rhoPythia( event[in2],
                      event[FinalPartPos[i]], event[in1], -1 );
      pt12 = min(pt12, temp);
    }

    if (withColour) {
      // Compute pythia FSR separation between two jets,
      // with knowledge of colour connections
      for(int j=0; j < int(FinalPartPos.size()); ++j) {

        // Find recoiler in event
        if ( i!=j ){
          bool isHard = false;
          int radAcl = event[FinalPartPos[i]].acol();
          int radCol = event[FinalPartPos[i]].col();
          int emtAcl = event[FinalPartPos[j]].acol();
          int emtCol = event[FinalPartPos[j]].col();
          int iRec = -1;
          // Check in final state
          if (iRec <= 0 && radAcl > 0 && radAcl != emtCol)
            iRec = findColour(radAcl, FinalPartPos[i], FinalPartPos[j],
                     event,1, isHard);
          if (iRec <= 0 && radCol > 0 && radCol != emtAcl)
            iRec = findColour(radCol, FinalPartPos[i], FinalPartPos[j],
                     event,1, isHard);
          if (iRec <= 0 && emtAcl > 0 && emtAcl != radCol)
            iRec = findColour(emtAcl, FinalPartPos[i], FinalPartPos[j],
                     event,1, isHard);
          if (iRec <= 0 && emtCol > 0 && emtCol != radAcl)
            iRec = findColour(emtCol, FinalPartPos[i], FinalPartPos[j],
                     event,1, isHard);

          // Check in initial state
          if (iRec <= 0 && radAcl > 0 && radAcl != emtCol)
            iRec = findColour(radAcl, FinalPartPos[i], FinalPartPos[j],
                     event,2, isHard);
          if (iRec <= 0 && radCol > 0 && radCol != emtAcl)
            iRec = findColour(radCol, FinalPartPos[i], FinalPartPos[j],
                     event,2, isHard);
          if (iRec <= 0 && emtAcl > 0 && emtAcl != radCol)
            iRec = findColour(emtAcl, FinalPartPos[i], FinalPartPos[j],
                     event,2, isHard);
          if (iRec <= 0 && emtCol > 0 && emtCol != radAcl)
            iRec = findColour(emtCol, FinalPartPos[i], FinalPartPos[j],
                     event,2, isHard);

          if (iRec > 0) {
            double temp = rhoPythia( event[FinalPartPos[i]],
                                     event[FinalPartPos[j]],
                                     event[iRec], 1 );
            pt12 = min(pt12, temp);
          }
        }
      }

    } else {
      // Compute pythia FSR separation between two jets,
      // without any knowledge of colour connections
      for(int j=0; j < int(FinalPartPos.size()); ++j) {
        for(int k=0; k < int(FinalPartPos.size()); ++k) {
          // Allow any parton as recoiler
          if ( (i != j) && (i != k) && (j != k) ){

            double temp = 0.;
            // Only check splittings allowed by flavour, e.g.
            // only q -> qg and g -> qqbar
            temp = rhoPythia( event[FinalPartPos[i]],
                              event[FinalPartPos[j]],
                              event[FinalPartPos[k]], 1 );
            pt12 = min(pt12, temp);
          }
        }
      }

      // Compute pythia FSR separation between two jets, with initial recoiler
      // without any knowledge of colour connections
      if ( event[in1].colType() != 0 && event[in2].colType() != 0) {
        for(int j=0; j < int(FinalPartPos.size()); ++j) {
          // Allow both initial partons as recoiler
          if ( i != j ){
            // Check with first initial as recoiler
            double temp = rhoPythia( event[FinalPartPos[i]],
                                     event[FinalPartPos[j]],
                                     event[in1], 1 );
            pt12 = min(pt12, temp);
            // Check with second initial as recoiler
            temp        = rhoPythia( event[FinalPartPos[i]],
                                     event[FinalPartPos[j]],
                                     event[in2], 1 );
            pt12 = min(pt12, temp);
          }
        }
      }

    }
    // Reset minimal y separation
    ptmin = min(ptmin,pt12);
  }

  return ptmin;

}

//--------------------------------------------------------------------------

// Function to compute "pythia pT separation" from Particle input, as a helper
// for rhoms(...).

double MergingHooks::rhoPythia(const Particle& RadAfterBranch,
              const Particle& EmtAfterBranch,
              const Particle& RecAfterBranch, int ShowerType){

  // Save type: 1 = FSR pT definition, else ISR definition
  int Type   = ShowerType;
  // Calculate virtuality of splitting
  int sign = (Type==1) ? 1 : -1;
  Vec4 Q(RadAfterBranch.p() + sign*EmtAfterBranch.p());
  double Qsq = sign * Q.m2Calc();
  // Mass term of radiator
  double m2Rad = ( includeMassive()
               && abs(RadAfterBranch.id()) >= 4
               && abs(RadAfterBranch.id()) < 7)
               ? pow(particleDataPtr->m0(RadAfterBranch.id()), 2)
               : 0.;
  // Construct 2->3 variables for FSR
  Vec4   sum     = RadAfterBranch.p() + RecAfterBranch.p()
                 + EmtAfterBranch.p();
  double m2Dip = sum.m2Calc();
  double x1 = 2. * (sum * RadAfterBranch.p()) / m2Dip;
  double x3 = 2. * (sum * EmtAfterBranch.p()) / m2Dip;
  // Construct momenta of dipole before/after splitting for ISR
  Vec4 qBR(RadAfterBranch.p() - EmtAfterBranch.p() + RecAfterBranch.p());
  Vec4 qAR(RadAfterBranch.p() + RecAfterBranch.p());
  // Calculate z of splitting, different for FSR and ISR
  double z = (Type==1) ? x1/(x1+x3)
                     : (qBR.m2Calc())/( qAR.m2Calc());
  // Separation of splitting, different for FSR and ISR
  double pTpyth = (Type==1) ? z*(1.-z) : (1.-z);
  // pT^2 = separation*virtuality
  pTpyth *= (Qsq - sign*m2Rad);
  if (pTpyth < 0.) pTpyth = 0.;
  // Return pT
  return sqrt(pTpyth);
}

//--------------------------------------------------------------------------

// Function to find a colour (anticolour) index in the input event.
// Helper for rhoms
// IN  int col       : Colour tag to be investigated
//     int iExclude1 : Identifier of first particle to be excluded
//                     from search
//     int iExclude2 : Identifier of second particle to be excluded
//                     from  search
//     Event event   : event to be searched for colour tag
//     int type      : Tag to define if col should be counted as
//                      colour (type = 1) [->find anti-colour index
//                                         contracted with col]
//                      anticolour (type = 2) [->find colour index
//                                         contracted with col]
// OUT int           : Position of particle in event record
//                     contraced with col [0 if col is free tag]

int MergingHooks::findColour(int col, int iExclude1, int iExclude2,
      const Event& event, int type, bool isHardIn){

  bool isHard = isHardIn;
  int index = 0;

  if (isHard){
    // Search event record for matching colour & anticolour
    for(int n = 0; n < event.size(); ++n) {
      if ( n != iExclude1 && n != iExclude2
        && event[n].colType() != 0
        &&(   event[n].status() > 0          // Check outgoing
           || event[n].status() == -21) ) {  // Check incoming
         if ( event[n].acol() == col ) {
          index = -n;
          break;
        }
        if ( event[n].col()  == col ){
          index =  n;
          break;
        }
      }
    }
  } else {

    // Search event record for matching colour & anticolour
    for(int n = 0; n < event.size(); ++n) {
      if (  n != iExclude1 && n != iExclude2
        && event[n].colType() != 0
        &&(   event[n].status() == 43        // Check outgoing from ISR
           || event[n].status() == 51        // Check outgoing from FSR
           || event[n].status() == 52        // Check outgoing from FSR
           || event[n].status() == -41       // first initial
           || event[n].status() == -42) ) {  // second initial
        if ( event[n].acol() == col ) {
          index = -n;
          break;
        }
        if ( event[n].col()  == col ){
          index =  n;
          break;
        }
      }
    }
  }
  // if no matching colour / anticolour has been found, return false
  if ( type == 1 && index < 0) return abs(index);
  if ( type == 2 && index > 0) return abs(index);

  return 0;
}

//--------------------------------------------------------------------------

// Find the if the event passes the Delta R_{ij}, pT_{i} and Q_{ij} cuts on
// the matrix element, as needed for cut-based merging scale definition.

double MergingHooks::cutbasedms( const Event& event ){

  // Only check first emission.
  if (!isFirstEmission(event)) return -1.;

  // Save allowed final state particles
  vector<int> partons;
  for( int i=0; i < event.size(); ++i) {
    if ( event[i].isFinal()
      && isInHard( i, event )
      && checkAgainstCut(event[i]) )
      partons.push_back(i);
  }

  // Declare overall veto
  bool doVeto = false;
  // Declare vetoes
  bool vetoPT  = false;
  bool vetoRjj = false;
  bool vetoMjj = false;
  // Declare cuts used in matrix element
  double pTjmin = pTiMS();
  double mjjmin = QijMS();
  double rjjmin = dRijMS();

  // Declare minimum values
  double minPT  = event[0].e();
  double minRJJ = 10.;
  double minMJJ = event[0].e();

  // Check matrix element cuts
  for( int i=0; i < int(partons.size()); ++i){
    // Save pT value
    minPT = min(minPT,event[partons[i]].pT());

    // Check two-parton cuts
    for( int j=0; j < int(partons.size()); ++j){
      if (i!=j){

        // Save delta R value
        minRJJ = min(minRJJ, deltaRij( event[partons[i]].p(),
              event[partons[j]].p()) );
        // Save delta R value
        minMJJ = min(minMJJ, ( event[partons[i]].p()
                              +event[partons[j]].p() ).mCalc() );

      }
    }
  // Done with cut evaluation
  }

  // Check if all partons are in the matrix element region
  if (minPT  > pTjmin) vetoPT  = true;
  if (minRJJ > rjjmin) vetoRjj = true;
  if (minMJJ > mjjmin) vetoMjj = true;

  // In the matrix element calculation, we have rejected the events if any of
  // the cuts had not been fulfilled,
  // i.e. we are in the matrix element domain if all cuts are fulfilled.
  // We veto any emission in the ME region.
  // Disregard the two-parton cuts if only one parton in the final state
  if (int(partons.size() == 1))
    doVeto = vetoPT;
  else
  // Veto if the combination of cuts would be in the ME region
    doVeto = vetoPT && vetoRjj && vetoMjj;

  // If event is above merging scale, veto
  if (doVeto) return 1.;

  // Else, do nothing
  return -1.;

}

//--------------------------------------------------------------------------

// Function to compute Delta R separation from 4-vector input.

double MergingHooks::deltaRij(Vec4 jet1, Vec4 jet2){

  // Declare return variable
  double deltaR = 0.;
  // Get delta_eta and cosh(Delta_eta) for hadronic collisions
  double eta1 = 0.5*log( (jet1.e() + jet1.pz()) / (jet1.e() - jet1.pz()) );
  double eta2 = 0.5*log( (jet2.e() + jet2.pz()) / (jet2.e() - jet2.pz()) );
  // Get delta_phi and cos(Delta_phi) for hadronic collisions
  double pt1 = sqrt( pow(jet1.px(),2) + pow(jet1.py(),2) );
  double pt2 = sqrt( pow(jet2.px(),2) + pow(jet2.py(),2) );
  double cosdPhi = ( jet1.px()*jet2.px() + jet1.py()*jet2.py() ) / (pt1*pt2);
  double dPhi = acos( cosdPhi );
  // Calculate kT durham like fastjet
  deltaR = sqrt(pow(eta1-eta2,2) + pow(dPhi,2));
  // Return kT
  return deltaR;
}

//==========================================================================

} // end namespace Pythia8
