// Settings.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the settings database.
// Flag: helper class with bool flags.
// Mode: helper class with int modes.
// Parm: (short for parameter) helper class with double parameters.
// Word: helper class with string words.
// MVec: vector of Modes (integers).
// PVec: vector of Parms (doubles).
// Settings: maps of flags, modes, parms and words with input/output.

#ifndef Pythia8_Settings_H
#define Pythia8_Settings_H

#include "Pythia8/Info.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Class for bool flags.

class Flag {

public:

  // Constructor
  Flag(string nameIn = " ", bool defaultIn = false) : name(nameIn),
    valNow(defaultIn) , valDefault(defaultIn) { }

  // Data members.
  string name;
  bool   valNow, valDefault;

};

//==========================================================================

// Class for integer modes.

class Mode {

public:

  // Constructor
  Mode(string nameIn = " ", int defaultIn = 0, bool hasMinIn = false,
    bool hasMaxIn = false, int minIn = 0,  int maxIn = 0,
    bool optOnlyIn = false) :  name(nameIn), valNow(defaultIn),
    valDefault(defaultIn), hasMin(hasMinIn), hasMax(hasMaxIn),
    valMin(minIn), valMax(maxIn), optOnly(optOnlyIn)  { }

  // Data members.
  string name;
  int    valNow, valDefault;
  bool   hasMin, hasMax;
  int    valMin, valMax;
  bool   optOnly;

};

//==========================================================================

// Class for double parms (where parm is shorthand for parameter).

class Parm {

public:

  // Constructor
  Parm(string nameIn = " ", double defaultIn = 0.,
    bool hasMinIn = false, bool hasMaxIn = false, double minIn = 0.,
    double maxIn = 0.) :  name(nameIn), valNow(defaultIn),
    valDefault(defaultIn), hasMin(hasMinIn), hasMax(hasMaxIn),
    valMin(minIn), valMax(maxIn) { }

  // Data members.
  string name;
  double valNow, valDefault;
  bool   hasMin, hasMax;
  double valMin, valMax;

};

//==========================================================================

// Class for string words.

class Word {

public:

  // Constructor
  Word(string nameIn = " ", string defaultIn = " ") : name(nameIn),
    valNow(defaultIn) , valDefault(defaultIn) { }

  // Data members.
  string name, valNow, valDefault;

};

//==========================================================================

// Class for vector of bool flags.

class FVec {

public:

  // Constructor
  FVec(string nameIn = " ", vector<bool> defaultIn = vector<bool>(1, false)) :
    name(nameIn), valNow(defaultIn) , valDefault(defaultIn) { }

  // Data members.
  string name;
  vector<bool> valNow, valDefault;

};

//==========================================================================

// Class for vector of integers.

class MVec {

public:

  // Constructor
  MVec(string nameIn = " ", vector<int> defaultIn = vector<int>(1, 0),
    bool hasMinIn = false, bool hasMaxIn = false, int minIn = 0,
    int maxIn = 0) :  name(nameIn), valNow(defaultIn),
    valDefault(defaultIn), hasMin(hasMinIn), hasMax(hasMaxIn),
    valMin(minIn), valMax(maxIn) { }

  // Data members.
  string name;
  vector<int> valNow, valDefault;
  bool   hasMin, hasMax;
  int    valMin, valMax;

};

//==========================================================================

// Class for vector of doubles.

class PVec {

public:

  // Constructor
  PVec(string nameIn = " ", vector<double> defaultIn = vector<double>(1, 0.),
    bool hasMinIn = false, bool hasMaxIn = false, double minIn = 0.,
    double maxIn = 0.) :  name(nameIn), valNow(defaultIn),
    valDefault(defaultIn), hasMin(hasMinIn), hasMax(hasMaxIn),
    valMin(minIn), valMax(maxIn) { }

  // Data members.
  string name;
  vector<double> valNow, valDefault;
  bool   hasMin, hasMax;
  double valMin, valMax;

};

//==========================================================================

// This class holds info on flags (bool), modes (int), parms (double),
// words (string), fvecs (vector of bool), mvecs (vector of int) and pvecs
// (vector of double).

class Settings {

public:

  // Constructor.
  Settings() : isInit(false), readingFailedSave(false) {}

  // Initialize Info pointer.
  void initPtr(Info* infoPtrIn) {infoPtr = infoPtrIn;}

  // Read in database from specific file.
  bool init(string startFile = "../xmldoc/Index.xml", bool append = false,
    ostream& os = cout) ;

  // Overwrite existing database by reading from specific file.
  bool reInit(string startFile = "../xmldoc/Index.xml", ostream& os = cout) ;

  // Read in one update from a single line.
  bool readString(string line, bool warn = true, ostream& os = cout) ;

  // Keep track whether any readings have failed, invalidating run setup.
  bool readingFailed() {return readingFailedSave;}

  // Write updates or everything to user-defined file.
  bool writeFile(string toFile, bool writeAll = false) ;
  bool writeFile(ostream& os = cout, bool writeAll = false) ;

  // Print out table of database, either all or only changed ones,
  // or ones containing a given string.
  void listAll(ostream& os = cout) {
    list( true, false, " ", os); }
  void listChanged(ostream& os = cout) {
    list (false, false, " ", os); }
  void list(string match, ostream& os = cout) {
    list (false, true, match, os); }

  // Reset all values to their defaults.
  void resetAll() ;

  // Query existence of an entry.
  bool isFlag(string keyIn) {
    return (flags.find(toLower(keyIn)) != flags.end()); }
  bool isMode(string keyIn) {
    return (modes.find(toLower(keyIn)) != modes.end()); }
  bool isParm(string keyIn) {
    return (parms.find(toLower(keyIn)) != parms.end()); }
  bool isWord(string keyIn) {
    return (words.find(toLower(keyIn)) != words.end()); }
  bool isFVec(string keyIn) {
    return (fvecs.find(toLower(keyIn)) != fvecs.end()); }
  bool isMVec(string keyIn) {
    return (mvecs.find(toLower(keyIn)) != mvecs.end()); }
  bool isPVec(string keyIn) {
    return (pvecs.find(toLower(keyIn)) != pvecs.end()); }

  // Add new entry.
  void addFlag(string keyIn, bool defaultIn) {
    flags[toLower(keyIn)] = Flag(keyIn, defaultIn); }
  void addMode(string keyIn, int defaultIn, bool hasMinIn,
    bool hasMaxIn, int minIn, int maxIn, bool optOnlyIn = false) {
    modes[toLower(keyIn)] = Mode(keyIn, defaultIn, hasMinIn, hasMaxIn,
    minIn, maxIn, optOnlyIn); }
  void addParm(string keyIn, double defaultIn, bool hasMinIn,
    bool hasMaxIn, double minIn, double maxIn) { parms[toLower(keyIn)]
    = Parm(keyIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn); }
  void addWord(string keyIn, string defaultIn) {
    words[toLower(keyIn)] = Word(keyIn, defaultIn); }
  void addFVec(string keyIn, vector<bool> defaultIn) {
    fvecs[toLower(keyIn)] = FVec(keyIn, defaultIn); }
  void addMVec(string keyIn, vector<int> defaultIn, bool hasMinIn,
    bool hasMaxIn, int minIn, int maxIn) { mvecs[toLower(keyIn)]
    = MVec(keyIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn); }
   void addPVec(string keyIn, vector<double> defaultIn, bool hasMinIn,
    bool hasMaxIn, double minIn, double maxIn) { pvecs[toLower(keyIn)]
    = PVec(keyIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn); }

  // Give back current value, with check that key exists.
  bool   flag(string keyIn);
  int    mode(string keyIn);
  double parm(string keyIn);
  string word(string keyIn);
  vector<bool>   fvec(string keyIn);
  vector<int>    mvec(string keyIn);
  vector<double> pvec(string keyIn);

  // Give back default value, with check that key exists.
  bool   flagDefault(string keyIn);
  int    modeDefault(string keyIn);
  double parmDefault(string keyIn);
  string wordDefault(string keyIn);
  vector<bool>   fvecDefault(string keyIn);
  vector<int>    mvecDefault(string keyIn);
  vector<double> pvecDefault(string keyIn);

  // Give back a map of all entries whose names match the string "match".
  map<string, Flag> getFlagMap(string match);
  map<string, Mode> getModeMap(string match);
  map<string, Parm> getParmMap(string match);
  map<string, Word> getWordMap(string match);
  map<string, FVec> getFVecMap(string match);
  map<string, MVec> getMVecMap(string match);
  map<string, PVec> getPVecMap(string match);

  // Change current value, respecting limits.
  void flag(string keyIn, bool nowIn);
  bool mode(string keyIn, int nowIn);
  void parm(string keyIn, double nowIn);
  void word(string keyIn, string nowIn);
  void fvec(string keyIn, vector<bool> nowIn);
  void mvec(string keyIn, vector<int> nowIn);
  void pvec(string keyIn, vector<double> nowIn);

  // Change current value, disregarding limits.
  void forceMode(string keyIn, int nowIn);
  void forceParm(string keyIn, double nowIn);
  void forceMVec(string keyIn, vector<int> nowIn);
  void forcePVec(string keyIn, vector<double> nowIn);

  // Restore current value to default.
  void resetFlag(string keyIn);
  void resetMode(string keyIn);
  void resetParm(string keyIn);
  void resetWord(string keyIn);
  void resetFVec(string keyIn);
  void resetMVec(string keyIn);
  void resetPVec(string keyIn);

private:

  // Pointer to various information on the generation.
  Info* infoPtr;

  // Map for bool flags.
  map<string, Flag> flags;

  // Map for integer modes.
  map<string, Mode> modes;

  // Map for double parms.
  map<string, Parm> parms;

  // Map for string words.
  map<string, Word> words;

  // Map for vectors of bool.
  map<string, FVec> fvecs;

  // Map for vectors of int.
  map<string, MVec> mvecs;

  // Map for vectors of double.
  map<string, PVec> pvecs;

  // Flags that initialization has been performed; whether any failures.
  bool isInit, readingFailedSave;

  // Print out table of database, called from listAll and listChanged.
  void list(bool doListAll, bool doListString, string match,
    ostream& os = cout);

  // Master switch for program printout.
  void printQuiet(bool quiet);

  // Restore settings used in tunes to e+e- and pp/ppbar data.
  void resetTuneEE();
  void resetTunePP();

  // Initialize tunes to e+e- and pp/ppbar data.
  void initTuneEE(int eeTune);
  void initTunePP(int ppTune);

  // Useful functions for string handling.
  string toLower(const string& name);
  bool   boolString(string tag);
  string attributeValue(string line, string attribute);
  bool   boolAttributeValue(string line, string attribute);
  int    intAttributeValue(string line, string attribute);
  double doubleAttributeValue(string line, string attribute);
  vector<bool>   boolVectorAttributeValue(string line, string attribute);
  vector<int>    intVectorAttributeValue(string line, string attribute);
  vector<double> doubleVectorAttributeValue(string line, string attribute);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Settings_H
