// SusyLesHouches.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for SUSY Les Houches Accord functionality
// This part of the SLHA interface basically contains the Pythia-independent
// SLHA read/write and processing utilities, which would be common to any
// SLHA interface.
// (The Pythia-specific components reside in the SLHAinterface class.)

#ifndef Pythia8_SLHA_H
#define Pythia8_SLHA_H

#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

//************************* SLHA AUX CLASSES *****************************//

  //class LHblock: the generic SLHA block (see below for matrices)
  //Explicit typing required, e.g. block<double> minpar;
  template <class T> class LHblock {

  public:

    //Constructor.
    LHblock<T>() : idnow(0), qDRbar(), i() {} ;

    //Does block exist?
    bool exists() { return int(entry.size()) == 0 ? false : true ; };
    //Clear block
    void clear() { entry.clear(); };

    //set: set block entry values.
    //Possible return values from set:
    // 0: normal return. Entry did not previously exist and has been created.
    // 1: normal return. Entry did previously exist and has been overwritten.
    //-1: failure.
    int set(int iIn,T valIn) {
      int alreadyexisting=exists(iIn)?1:0;
      entry[iIn]=valIn;
      return alreadyexisting;
    };
    // Read index and value from SLHA data line
    int set(istringstream& linestream, bool indexed=true) {
      i = 0;
      if (indexed) linestream >> i >> val;
      else linestream >> val;
      return linestream ? set(i,val) : -1;
    };
    // With i already given, read value from remaining SLHA data line
    int set(int iIn,istringstream& linestream) {
      linestream >> val;
      return linestream ? set(iIn,val) : -1;
    };
    // Shorthand for entry[0]. Used e.g. for block ALPHA.
    void set(T valIn) { entry[0]=valIn; };

    // Does entry i already exist in this block?
    bool exists(int iIn) {return entry.find(iIn) != entry.end()
      ? true : false;};

    // Indexing with (). Output only.
    T operator()() {
      if (exists(0)) {return entry[0];} else {T dummy(0); return dummy;};
    };
    T operator()(int iIn) {
      if (exists(iIn)) {return entry[iIn];} else {T dummy(0); return dummy;};
    };

    // Size of map
    int size() {return int(entry.size());};

    // First and next key code
    int first() { idnow = entry.begin()->first; return idnow; };
    int next() {
      typename map<int,T>::iterator itnow;
      itnow = ++entry.find(idnow);
      if ( itnow == entry.end() ) itnow=entry.begin();
      return idnow = itnow->first;
    };

    // Simple print utility
    void list() {
      bool finished=false;
      int ibegin=first();
      i=ibegin;
      while (!finished) {
        cout << "  "<< i << " " << entry[i] <<endl;
        i=next();
        if (i == ibegin) finished=true;
      };
    };

    // Special for DRbar running blocks.
    void setq(double qIn) { qDRbar=qIn; }
    double q() { return qDRbar; }

  protected:
    map<int,T> entry;

  private:
    int idnow;
    double qDRbar;
    //Auxiliary vars
    int i;
    T val;
  };

  // Derived class for generic blocks containing vectors of strings.
  class LHgenericBlock : public LHblock<string> {

  public:

    //Constructor.
    LHgenericBlock() { } ;

    // Read index and value from SLHA data line
    int set(string lineIn) {
      entry[entry.size()] = lineIn;
      return 0;
    };

  };

  // class LHmatrixBlock: the generic SLHA matrix
  // Explicit sizing required, e.g.LHmatrixBlock<4> nmix;
  template <int size> class LHmatrixBlock {
  public:
    //Constructor. Set uninitialized and explicitly zero.
    LHmatrixBlock<size>() : entry(), qDRbar(), val() {
      initialized=false;
      for (i=1;i<=size;i++) {
        for (j=1;j<=size;j++) {
          entry[i][j]=0.0;
        };
      };
    };

    // Assignment
    LHmatrixBlock& operator=(const LHmatrixBlock& m) {
      if (this != &m) {
        for (i=0;i<size;i++) for (j=0;j<=size;j++) entry[i][j] = m(i,j);
        qDRbar = m.qDRbar;
        initialized = m.initialized;
      }
      return *this; };

    // Does this matrix contain any entries?
    bool exists() { return initialized; };
    // Clear initialized flag
    void clear() { initialized=false; };

    // Set matrix entry
    int set(int iIn,int jIn, double valIn) {
      if (iIn>0 && jIn>0 && iIn<=size && jIn<=size) {
        entry[iIn][jIn]=valIn;
        initialized=true;
        return 0;
      } else {
        return -1;
      };
    };

    // Set entry from linestream (used during file read)
    int set(istringstream& linestream) {
      linestream >> i >> j >> val;
      return linestream ? set(i,j,val) : -1;
    };

    // () Overloading: Get entry
    double operator()(int iIn, int jIn) const {
      return (iIn <= size && jIn <= size && iIn > 0 && jIn > 0) ?
        entry[iIn][jIn] : 0.0;
    };

    // Set and get scale for DRbar running LHblocks.
    void setq(double qIn) { qDRbar=qIn; }
    double q() { return qDRbar; }

    // Simple print utility, to be elaborated on.
    void list() {
      for (i=1;i<=size;i++) {
        cout << "   "<<i << " " ;
        for (j=1;j<=size;j++) cout << entry[i][j] << " ";
        cout << endl;
      };
    };

  private:
    bool initialized;
    double entry[size+1][size+1];
    double qDRbar;
    //Auxiliary vars
    int i,j;
    double val;
  };

  // class tensorBlock: the generic SLHA tensor
  // Explicit sizing required, e.g. tensorBlock<3> rvlam;
  template <int size> class LHtensor3Block {
  public:
    //Constructor. Set uninitialized and explicitly zero.
    LHtensor3Block<size>() : entry(), qDRbar(), val() {
      initialized=false;
      for (i=1;i<=size;i++) {
        for (j=1;j<=size;j++) {
          for (k=1;k<=size;k++) {
            entry[i][j][k]=0.0;
          };
        };
      };
    };

    // Assignment
    LHtensor3Block& operator=(const LHtensor3Block& m) {
      if (this != &m) {
        for (i=0;i<size;i++) for (j=0;j<=size;j++) for (k=0;k<=size;k++)
          entry[i][j][k] = m(i,j,k);
        qDRbar = m.qDRbar;
        initialized = m.initialized;
      }
      return *this; };

    // Does this matrix contain any entries?
    bool exists() { return initialized; };
    // Clear initialized flag
    void clear() { initialized=false; };

    // Set matrix entry
    int set(int iIn,int jIn, int kIn, double valIn) {
      if (iIn>0 && jIn>0 && kIn>0 && iIn<=size && jIn<=size && kIn<=size) {
        entry[iIn][jIn][kIn]=valIn;
        initialized=true;
        return 0;
      } else {
        return -1;
      };
    };

    // Set entry from linestream (used during file read)
    int set(istringstream& linestream) {
      linestream >> i >> j >> k >> val;
      return linestream ? set(i,j,k,val) : -1;
    };

    // () Overloading: Get entry
    double operator()(int iIn, int jIn, int kIn) const {
      return (iIn <= size && jIn <= size && kIn <= size && iIn > 0
        && jIn > 0 && kIn > 0) ? entry[iIn][jIn][kIn] : 0.0;
    };

    // Set and get scale for DRbar running LHblocks.
    void setq(double qIn) { qDRbar=qIn; }
    double q() { return qDRbar; }

    // Simple print utility, to be elaborated on.
    void list() {
      for (i=1;i<=size;i++) {
        for (j=1;j<=size;j++) {
          cout << "   "<<i << " "<<j << " " ;
          for (k=1;k<=size;k++) {
            cout << entry[i][j][k] << " ";
            cout << endl;
          };
        };
      };
    };

  private:
    bool initialized;
    double entry[size+1][size+1][size+1];
    double qDRbar;
    //Auxiliary vars
    int i,j,k;
    double val;
  };

  //*************************** DECAY TABLES ***************************//

  class LHdecayChannel {
  public:

    LHdecayChannel() : brat(0.0) {};
    LHdecayChannel(double bratIn, int nDaIn, vector<int> idDaIn,
      string cIn="") : brat() { setChannel(bratIn,nDaIn,idDaIn,cIn);
    }

    // Functions to set decay channel information
    void setChannel(double bratIn, int nDaIn, vector<int> idDaIn,
      string cIn="") {
      brat    = bratIn;
      for (int i=0; i<=nDaIn; i++) {
        if (i < int(idDaIn.size())) idDa.push_back(idDaIn[i]);
        comment = cIn;
      }
    }
    void setBrat(double bratIn) {brat=bratIn;}
    void setIdDa(vector<int> idDaIn) {idDa = idDaIn;}

    // Functions to get decay channel information
    double getBrat() {return brat;}
    int getNDa() {return int(idDa.size());}
    vector<int> getIdDa() {return idDa;}
    string getComment() {return comment;}

  private:
    double brat;
    vector<int> idDa;
    string comment;

  };

  class LHdecayTable {
  public:

  LHdecayTable() : id(0), width(0.0) {};
  LHdecayTable(int idIn) : id(idIn), width(0.0) {};
  LHdecayTable(int idIn, double widthIn) : id(idIn), width(widthIn) {};

    // Functions to get PDG code (id) and width
    int    getId() {return id;}
    double getWidth() {return width;}

    // Functions to set PDG code (id) and width
    void setId(int idIn) {id = idIn;}
    void setWidth(double widthIn) {width=widthIn;}

    // Function to reset size and width (width -> 0 by default)
    void reset(double widthIn=0.0) {table.resize(0); width=widthIn;}

    // Function to add another decay channel
    void addChannel(LHdecayChannel channelIn) {table.push_back(channelIn);}
    void addChannel(double bratIn, int nDaIn, vector<int> idDaIn,
      string cIn="") {
      LHdecayChannel newChannel(bratIn, nDaIn, idDaIn, cIn);
      table.push_back(newChannel);
    }

    // Function to return number of decay channels
    int size() {return int(table.size());}

    // Function to return a branching ratio
    double getBrat(int iChannel) {
      if (iChannel >= 0 && iChannel < int(table.size())) {
        return table[iChannel].getBrat();
      } else {
        return 0.0;
      }
    }
    // Function to return daughter PDG codes
    vector<int> getIdDa(int iChannel) {
      if (iChannel >= 0 && iChannel < int(table.size())) {
        return table[iChannel].getIdDa();
      } else {
        vector<int> dum;
        return dum;
      }
    }
    // Function to return a decay channel
    LHdecayChannel getChannel(int iChannel) {
      if (iChannel >= 0 && iChannel < int(table.size())) {
        return table[iChannel];
      } else {
        LHdecayChannel dum;
        return dum;
      }
    }

  private:
    int id;
    double width;
    vector<LHdecayChannel> table;

  };

//==========================================================================

class SusyLesHouches {

public:

  //Constructor, with and without filename.
  SusyLesHouches(int verboseIn=1) : verboseSav(verboseIn),
    headerPrinted(false), footerPrinted(false), filePrinted(false),
    slhaRead(false), lhefRead(false), lhefSlha(false), useDecay(true) {};
  SusyLesHouches(string filename, int verboseIn=1) : verboseSav(verboseIn),
    headerPrinted(false), footerPrinted(false), filePrinted(false),
    slhaRead(true), lhefRead(false), lhefSlha(false), useDecay(true)
    {readFile(filename);};

  //***************************** SLHA FILE I/O *****************************//
  // Read and write SLHA files
  int readFile(string slhaFileIn="slha.spc",int verboseIn=1,
    bool useDecayIn=true);
  int readFile(istream& ,int verboseIn=1,
    bool useDecayIn=true);
  //int writeFile(string filename): write SLHA file on filename

  //Output utilities
  void listHeader();   // print Header
  void listFooter();   // print Footer
  void listSpectrum(int ifail=0); // print Spectrum

  // Check spectrum and decays
  int checkSpectrum();

  // File Name (can be either SLHA or LHEF)
  string slhaFile;

  // Class for SLHA data entry
  class Entry {

  public:
    //Constructor.
    Entry() : isIntP(false), isDoubleP(false),
      isStringP(false), n(0), d(0.0), s(""), commentP("") {}

    // Generic functions to inquire whether an int, double, or string
    bool isInt(){return isIntP;}
    bool isDouble(){return isDoubleP;}
    bool isString(){return isStringP;}

    // = Overloading: Set entry to int, double, or string
    Entry& operator=(double& val)  {
      d=val;isIntP=false;isDoubleP=true;isStringP=false;
      return *this;
    };
    Entry& operator=(int& val)  {
      n=val;isIntP=true;isDoubleP=false;isStringP=false;
      return *this;
    };
    Entry& operator=(string& val)  {
      s=val;isIntP=false;isDoubleP=false;isStringP=true;
      return *this;
    };

    // Set and Get comment
    void setComment(string comment) {commentP=comment;}
    void getComment(string comment) {comment=commentP;}

    // Generic functions to get value
    bool get(int& val) {val=n; return isIntP;}
    bool get(double& val) {val=d; return isDoubleP;}
    bool get(string& val) {val=s; return isStringP;}

  private:
    bool isIntP, isDoubleP, isStringP;
    int n;
    double d;
    string s;
    string commentP;

  };

  //*************************** THE SLHA1 BLOCKS ***************************//
  //Blocks for model definition:
  LHblock<int> modsel;
  LHblock<int> modsel21;
  LHblock<double> modsel12;
  LHblock<double> minpar;
  LHblock<double> extpar;
  LHblock<double> sminputs;
  //Blocks for RGE program specific output
  LHblock<string> spinfo;
  LHblock<string> spinfo3;
  LHblock<string> spinfo4;
  //Blocks for DCY program specific output
  LHblock<string> dcinfo;
  LHblock<string> dcinfo3;
  LHblock<string> dcinfo4;
  //Blocks for mass and coupling spectrum
  LHblock<double> mass;
  LHmatrixBlock<4> nmix;
  LHmatrixBlock<2> umix;
  LHmatrixBlock<2> vmix;
  LHmatrixBlock<2> stopmix;
  LHmatrixBlock<2> sbotmix;
  LHmatrixBlock<2> staumix;
  LHblock<double> alpha;
  LHblock<double> hmix;
  LHblock<double> gauge;
  LHblock<double> msoft;
  LHmatrixBlock<3> au;
  LHmatrixBlock<3> ad;
  LHmatrixBlock<3> ae;
  LHmatrixBlock<3> yu;
  LHmatrixBlock<3> yd;
  LHmatrixBlock<3> ye;

  //************************ THE SLHA1 DECAY TABLES ************************//
  vector<LHdecayTable> decays;
  map<int,int> decayIndices;

  //********************* THE BSM-SLHA QNUMBERS BLOCKS *********************//
  vector< LHblock<double> > qnumbers;     // Zero'th entry is PDG code
  vector< string > qnumbersName;
  vector< string > qnumbersAntiName;

  //*************************** THE SLHA2 BLOCKS ***************************//
  //Additions to SLHA1
  LHblock<double> qextpar;

  //FLV Input
  LHblock<double> vckmin;  // The input CKM Wolfenstein parms.
  LHblock<double> upmnsin; // The input PMNS PDG parms.
  LHmatrixBlock<3> msq2in; // The input upper off-diagonal msq2
  LHmatrixBlock<3> msu2in; // The input upper off-diagonal msu2
  LHmatrixBlock<3> msd2in; // The input upper off-diagonal msd2
  LHmatrixBlock<3> msl2in; // The input upper off-diagonal msl2
  LHmatrixBlock<3> mse2in; // The input upper off-diagonal mse2
  LHmatrixBlock<3> tuin;   // The input upper off-diagonal TU
  LHmatrixBlock<3> tdin;   // The input upper off-diagonal TD
  LHmatrixBlock<3> tein;   // The input upper off-diagonal TE
  //FLV Output
  LHmatrixBlock<3> vckm;    // The output DRbar running Re{VCKM} at Q
  LHmatrixBlock<3> upmns;   // The output DRbar running Re{UPMNS} at Q
  LHmatrixBlock<3> msq2;    // The output DRbar running msq2 at Q
  LHmatrixBlock<3> msu2;    // The output DRbar running msu2 at Q
  LHmatrixBlock<3> msd2;    // The output DRbar running msd2 at Q
  LHmatrixBlock<3> msl2;    // The output DRbar running msl2 at Q
  LHmatrixBlock<3> mse2;    // The output DRbar running mse2 at Q
  LHmatrixBlock<3> tu;      // The output DRbar running TU at Q
  LHmatrixBlock<3> td;      // The output DRbar running TD at Q
  LHmatrixBlock<3> te;      // The output DRbar running TE at Q
  LHmatrixBlock<6> usqmix;  // The Re{} up squark mixing matrix
  LHmatrixBlock<6> dsqmix;   // The Re{} down squark mixing matrix
  LHmatrixBlock<6> selmix;   // The Re{} selectron mixing matrix
  LHmatrixBlock<3> snumix;   // The Re{} sneutrino mixing matrix
  LHmatrixBlock<3> snsmix;   // The scalar sneutrino mixing matrix
  LHmatrixBlock<3> snamix;   // The pseudoscalar neutrino mixing matrix

  //RPV Input
  LHtensor3Block<3> rvlamllein; // The input LNV lambda couplings
  LHtensor3Block<3> rvlamlqdin; // The input LNV lambda' couplings
  LHtensor3Block<3> rvlamuddin; // The input BNV lambda'' couplings
  LHtensor3Block<3> rvtllein;   // The input LNV T couplings
  LHtensor3Block<3> rvtlqdin;   // The input LNV T' couplings
  LHtensor3Block<3> rvtuddin;   // The input BNV T'' couplings
  LHblock<double> rvkappain;    // The input LNV kappa couplings
  LHblock<double> rvdin;        // The input LNV D terms
  LHblock<double> rvm2lh1in;    // The input LNV m2LH1 couplings
  LHblock<double> rvsnvevin;    // The input LNV sneutrino vevs
  //RPV Output
  LHtensor3Block<3> rvlamlle;   // The output LNV lambda couplings
  LHtensor3Block<3> rvlamlqd;   // The output LNV lambda' couplings
  LHtensor3Block<3> rvlamudd;   // The output BNV lambda'' couplings
  LHtensor3Block<3> rvtlle;     // The output LNV T couplings
  LHtensor3Block<3> rvtlqd;     // The output LNV T' couplings
  LHtensor3Block<3> rvtudd;     // The output BNV T'' couplings
  LHblock<double> rvkappa;      // The output LNV kappa couplings
  LHblock<double> rvd;          // The output LNV D terms
  LHblock<double> rvm2lh1;      // The output LNV m2LH1 couplings
  LHblock<double> rvsnvev;      // The output LNV sneutrino vevs
  LHmatrixBlock<7> rvnmix;      // The RPV neutralino mixing matrix
  LHmatrixBlock<5> rvumix;      // The RPV chargino L mixing matrix
  LHmatrixBlock<5> rvvmix;      // The RPV chargino R mixing matrix
  LHmatrixBlock<5> rvhmix;      // The RPV neutral scalar mixing matrix
  LHmatrixBlock<5> rvamix;      // The RPV neutral pseudoscalar mixing matrix
  LHmatrixBlock<8> rvlmix;      // The RPV charged fermion mixing matrix

  //CPV Input
  LHblock<double> imminpar;
  LHblock<double> imextpar;
  //CPV Output
  LHmatrixBlock<4> cvhmix;   // The CPV Higgs mixing matrix
  LHmatrixBlock<4> imcvhmix; // Optional: imaginary components
  LHmatrixBlock<3> imau,imad,imae; // Im{} of AU, AD, AE
  LHblock<double> imhmix;
  LHblock<double> immsoft;

  //CPV + FLV Input
  LHmatrixBlock<3> immsq2in;  // The Im{} input upper off-diagonal msq2
  LHmatrixBlock<3> immsu2in;  // The Im{} input upper off-diagonal msu2
  LHmatrixBlock<3> immsd2in;  // The Im{} input upper off-diagonal msd2
  LHmatrixBlock<3> immsl2in;  // The Im{} input upper off-diagonal msl2
  LHmatrixBlock<3> immse2in;  // The Im{} input upper off-diagonal mse2
  LHmatrixBlock<3> imtuin,imtdin,imtein; // The Im{} input upper off-diagonal T
  //CPV + FLV Output
  LHmatrixBlock<3> imvckm;  // The output DRbar running Im{VCKM} at Q
  LHmatrixBlock<3> imupmns; // The output DRbar running Im{UPMNS} at Q
  LHmatrixBlock<3> immsq2;  // The output DRbar running msq2 at Q
  LHmatrixBlock<3> immsu2;  // The output DRbar running msu2 at Q
  LHmatrixBlock<3> immsd2;  // The output DRbar running msd2 at Q
  LHmatrixBlock<3> immsl2;  // The output DRbar running msl2 at Q
  LHmatrixBlock<3> immse2;  // The output DRbar running mse2 at Q
  LHmatrixBlock<3> imtu,imtd,imte; // Im{} of TU, TD, TE
  LHmatrixBlock<6> imusqmix;// The Im{} up squark mixing matrix
  LHmatrixBlock<6> imdsqmix; // The Im{} down squark mixing matrix
  LHmatrixBlock<6> imselmix; // The Im{} selectron mixing matrix
  LHmatrixBlock<3> imsnumix; // The Im{} sneutrino mixing matrix
  LHmatrixBlock<4> imnmix;   // The Im{} neutralino mixing matrix
  LHmatrixBlock<4> imumix;   // The Im{} chargino L mixing matrix
  LHmatrixBlock<4> imvmix;   // The Im{} chargino R mixing matrix

  //NMSSM Input
  //    All input is in EXTPAR
  //NMSSM Output
  LHblock<double> nmssmrun;  // The LHblock of NMSSM running parameters
  LHmatrixBlock<3> nmhmix;   // The NMSSM scalar Higgs mixing
  LHmatrixBlock<3> nmamix;   // The NMSSM pseudoscalar Higgs mixing
  LHmatrixBlock<5> nmnmix;   // The NMSSM neutralino mixing
  LHmatrixBlock<5> imnmnmix; //   Im{} (for future use)

  //*************************** SET BLOCK VALUE ****************************//
  template <class T> int set(string,T);
  template <class T> int set(string,int,T);
  template <class T> int set(string,int,int,T);
  template <class T> int set(string,int,int,int,T);

  //********************* GENERIC/USER-DEFINED BLOCKS **********************//
  // bool getEntry(name, indices, value)
  //      = true if LHblock and entry exists (value returned in value,
  //        typecast by user in call)
  //      = false otherwise
  map<string, LHgenericBlock> genericBlocks;
  template <class T> bool getEntry(string, T&);
  template <class T> bool getEntry(string, int, T&);
  template <class T> bool getEntry(string, int, int, T&);
  template <class T> bool getEntry(string, int, int, int, T&);
  template <class T> bool getEntry(string, vector<int>, T&);

  // Access/change verbose setting
  int verbose() {return verboseSav;}
  void verbose(int verboseIn) {verboseSav = verboseIn;}

  // Output of messages from SLHA interface
  void message(int, string,string ,int line=0);

  //***************************** SLHA PRIVATE *****************************//
private:
  //SLHA I/O
  int verboseSav;
  bool headerPrinted, footerPrinted, filePrinted;
  bool slhaRead, lhefRead, lhefSlha, useDecay;

};

//--------------------------------------------------------------------------

// utilities to set generic blocks

template <class T> int SusyLesHouches::set(string blockName, T val) {

  // Make sure everything is interpreted as lower case (for safety)
  toLowerRep(blockName);

  // Add new generic block if not already existing
  if (genericBlocks.find(blockName) == genericBlocks.end()) {
    LHgenericBlock gBlock;
    genericBlocks[blockName]=gBlock;
  }

  // Convert input value to string
  ostringstream lineStream;
  lineStream << val;
  return genericBlocks[blockName].set(lineStream.str());

}

template <class T> int SusyLesHouches::set(string blockName, int indx, T val) {

  // Make sure everything is interpreted as lower case (for safety)
  toLowerRep(blockName);

  // Add new generic block if not already existing
  if (genericBlocks.find(blockName) == genericBlocks.end()) {
    LHgenericBlock gBlock;
    genericBlocks[blockName]=gBlock;
  }

  // Convert input value to string
  ostringstream lineStream;
  lineStream << indx<<" "<<val;
  return genericBlocks[blockName].set(lineStream.str());

}

template <class T> int SusyLesHouches::set(string blockName, int indx,
                                           int jndx, T val) {

  // Make sure everything is interpreted as lower case (for safety)
  toLowerRep(blockName);

  // Add new generic block if not already existing
  if (genericBlocks.find(blockName) == genericBlocks.end()) {
    LHgenericBlock gBlock;
    genericBlocks[blockName]=gBlock;
  }

  // Convert input value to string
  ostringstream lineStream;
  lineStream << indx<<" "<<jndx<<" "<<val;
  return genericBlocks[blockName].set(lineStream.str());

}

template <class T> int SusyLesHouches::set(string blockName, int indx,
                                           int jndx, int kndx, T val) {

  // Make sure everything is interpreted as lower case (for safety)
  toLowerRep(blockName);

  // Add new generic block if not already existing
  if (genericBlocks.find(blockName) == genericBlocks.end()) {
    LHgenericBlock gBlock;
    genericBlocks[blockName]=gBlock;
  }

  // Convert input value to string
  ostringstream lineStream;
  lineStream << indx<<" "<<jndx<<" "<<kndx<<" "<<val;
  return genericBlocks[blockName].set(lineStream.str());

}

// utilities to read generic blocks

template <class T> bool SusyLesHouches::getEntry(string blockName, T& val) {

  // Make sure everything is interpret as lower case (for safety)
  toLowerRep(blockName);

  // Safety checks
  if (genericBlocks.find(blockName) == genericBlocks.end()) {
    message(1,"getEntry","attempting to extract entry from non-existent block "
            +blockName);
    return false;
  }
  if (genericBlocks[blockName].size() == 0) {
    message(1,"getEntry","attempting to extract entry from zero-size block "
            +blockName);
    return false;
  }
  if (genericBlocks[blockName].size() >= 2) {
    message(1,"getEntry","attempting to extract un-indexed entry "
      "from multi-entry block "+blockName);
    return false;
  }
  // Attempt to extract value as class T
  LHgenericBlock block = genericBlocks[blockName];
  istringstream linestream(block(0));
  linestream >> val;
  if ( !linestream ) {
    message(1,"getEntry","problem extracting un-indexed entry "
      "from block "+blockName);
    return false;
  }
  // If made it all the way here, value was successfully extracted.
  // Return true.
  return true;
}

template <class T> bool SusyLesHouches::getEntry(string blockName, int indx,
                                                 T& val) {

  // Make sure everything is interpret as lower case (for safety)
  toLowerRep(blockName);

  // Safety checks
  if (genericBlocks.find(blockName) == genericBlocks.end()) {
    message(1,"getEntry","attempting to extract entry from non-existent block "
            +blockName);
    return false;
  }
  if (genericBlocks[blockName].size() == 0) {
    message(1,"getEntry","attempting to extract entry from zero-size block "
            +blockName);
    return false;
  }
  // Attempt to extract indexed value as class T
  LHgenericBlock block = genericBlocks[blockName];
  // Loop over block contents, search for indexed entry with index i
  for (int jEntry = 0; jEntry < block.size(); jEntry++) {
    istringstream linestream(block(jEntry));
    // Buffer line according to format selected by T
    int indxNow;
    T valNow;
    linestream >> indxNow >> valNow;
    // If index found and value was readable, return true
    if (linestream && indxNow == indx) {
      val = valNow;
      return true;
    }
  }
  // If index not found or unreadable, return false
  message(1,"getEntry","problem extracting indexed entry from block "
          +blockName);
  return false;
}

template <class T> bool SusyLesHouches::getEntry(string blockName, int indx,
                                                 int jndx, T& val) {

  // Make sure everything is interpret as lower case (for safety)
  toLowerRep(blockName);

  // Safety checks
  if (genericBlocks.find(blockName) == genericBlocks.end()) {
    message(1,"getEntry","attempting to extract entry from non-existent block "
            +blockName);
    return false;
  }
  if (genericBlocks[blockName].size() == 0) {
    message(1,"getEntry","attempting to extract entry from zero-size block "
            +blockName);
    return false;
  }
  // Attempt to extract matrix-indexed value as class T
  LHgenericBlock block = genericBlocks[blockName];
  // Loop over block contents, search for indexed entry with indices i, j
  for (int jEntry = 0; jEntry < block.size(); jEntry++) {
    istringstream linestream(block(jEntry));
    // Buffer line according to format selected by T
    int indxNow, jndxNow;
    T valNow;
    linestream >> indxNow >> jndxNow >> valNow;
    // If index found and value was readable, return true
    if (linestream && indxNow == indx && jndxNow == jndx) {
      val = valNow;
      return true;
    }
  }
  // If index not found or unreadable, return false
  message(1,"getEntry","problem extracting matrix-indexed entry from block "
          +blockName);
  return false;
}

template <class T> bool SusyLesHouches::getEntry(string blockName, int indx,
                                                 int jndx, int kndx, T& val) {

  // Make sure everything is interpret as lower case (for safety)
  toLowerRep(blockName);

  // Safety checks
  if (genericBlocks.find(blockName) == genericBlocks.end()) {
    message(1,"getEntry","attempting to extract entry from non-existent block "
            +blockName);
    return false;
  }
  if (genericBlocks[blockName].size() == 0) {
    message(1,"getEntry","attempting to extract entry from zero-size block "
            +blockName);
    return false;
  }
  // Attempt to extract tensor-indexed value as class T
  LHgenericBlock block = genericBlocks[blockName];
  // Loop over block contents, search for indexed entry with indices i, j, k
  for (int jEntry = 0; jEntry < block.size(); jEntry++) {
    istringstream linestream(block(jEntry));
    // Buffer line according to format selected by T
    int indxNow, jndxNow, kndxNow;
    T valNow;
    linestream >> indxNow >> jndxNow >> kndxNow >> valNow;
    // If index found and value was readable, return true
    if (linestream && indxNow == indx && jndxNow == jndx && kndxNow == kndx) {
      val = valNow;
      return true;
    }
  }
  // If index not found or unreadable, return false
  message(1,"getEntry","problem extracting tensor-indexed entry from block "
          +blockName);
  return false;
 }

//==========================================================================


} // end of namespace Pythia8

#endif // end Pythia8_SLHA_H
