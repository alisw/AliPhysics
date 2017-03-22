// Settings.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Settings class.

#include "Pythia8/Settings.h"

// Allow string and character manipulation.
#include <cctype>

namespace Pythia8 {

//==========================================================================

// Settings class.
// This class contains flags, modes, parms and words used in generation.

//--------------------------------------------------------------------------

// Read in database from specific file.

bool Settings::init(string startFile, bool append, ostream& os) {

  // Don't initialize if it has already been done and not in append mode.
  if (isInit && !append) return true;
  int nError = 0;

  // List of files to be checked. Start with input file.
  vector<string> files;
  files.push_back(startFile);

  // If nontrivial startfile path, then use that for other files as well.
  string pathName = "";
  if (startFile.rfind("/") != string::npos)
    pathName = startFile.substr(0, startFile.rfind("/") + 1);

  // Loop over files. Open them for read.
  for (int i = 0; i < int(files.size()); ++i) {
    const char* cstring = files[i].c_str();
    ifstream is(cstring);

    // Check that instream is OK.
    if (!is.good()) {
      os << "\n PYTHIA Error: settings file " << files[i]
         << " not found" << endl;
      return false;
    }

    // Read in one line at a time.
    string line;
    while ( getline(is, line) ) {

      // Get first word of a line, to interpret it as tag.
      istringstream getfirst(line);
      string tag;
      getfirst >> tag;

      // Skip ahead if not interesting. Only look for new files in startfile.
      if (tag != "<flag" && tag != "<flagfix" && tag != "<mode"
         && tag != "<modeopen" && tag != "<modepick" && tag != "<modefix"
         && tag != "<parm" && tag != "<parmfix" && tag != "<word"
         && tag != "<wordfix" && tag != "<fvec" && tag != "<fvecfix"
         && tag != "<mvec" && tag != "<mvecfix"
         && tag != "<pvec" && tag != "<pvecfix" && tag != "<aidx") continue;

      // Read and append continuation line(s) if line does not contain >.
      while (line.find(">") == string::npos) {
        string addLine;
        getline(is, addLine);
        line += " " + addLine;
      }

      // Remove extra blanks before an = sign.
      while (line.find(" =") != string::npos) line.erase( line.find(" ="), 1);

      // Add file also to be read.
      if (tag == "<aidx") {
        string name = attributeValue( line, "href");
        if (name == "") {
          os << " PYTHIA Error: failed to find name attribute in line "
             << line << endl;
          ++nError;
          continue;
        }
        files.push_back(pathName + name + ".xml");
        continue;
      }

      // Find name attribute.
      string name = attributeValue( line, "name=");
      if (name == "") {
        os << " PYTHIA Error: failed to find name attribute in line "
           << line << endl;
        ++nError;
        continue;
      }

      // Check that default value attribute present, and whether max and min.
      if (line.find("default=") == string::npos) {
        os << " PYTHIA Error: failed to find default value token in line "
           << line << endl;
        ++nError;
        continue;
      }
      bool hasMin = (line.find("min=") != string::npos);
      bool hasMax = (line.find("max=") != string::npos);

      // Check for occurence of a bool and add to flag map.
      if (tag == "<flag" || tag == "<flagfix") {
        bool value = boolAttributeValue( line, "default=");
        addFlag( name, value);

      // Check for occurence of an int and add to mode map.
      } else if (tag == "<mode" || tag == "<modeopen"
        || tag == "<modepick" || tag == "<modefix") {
        int value    = intAttributeValue( line, "default=");
        int minVal   = intAttributeValue( line, "min=");
        int maxVal   = intAttributeValue( line, "max=");
        // Enforce check that only allowed options are accepted.
        bool optOnly = false;
        if (tag == "<modepick" && hasMin && hasMax) optOnly = true;
        if (tag == "<modefix") {
          hasMin  = true;
          hasMax  = true;
          minVal  = value;
          maxVal  = value;
          optOnly = true;
        }
        addMode( name, value, hasMin, hasMax, minVal, maxVal, optOnly);

      // Check for occurence of a double and add to parm map.
      } else if (tag == "<parm" || tag == "<parmfix") {
        double value  = doubleAttributeValue( line, "default=");
        double minVal = doubleAttributeValue( line, "min=");
        double maxVal = doubleAttributeValue( line, "max=");
        addParm( name, value, hasMin, hasMax, minVal, maxVal);

      // Check for occurence of a string and add to word map.
      } else if (tag == "<word" || tag == "<wordfix") {
        string value = attributeValue( line, "default=");
        addWord( name, value);

      // Check for occurence of a bool vector and add to fvec map.
      } else if (tag == "<fvec" || tag == "<fvecfix") {
        vector<bool> value = boolVectorAttributeValue( line, "default=");
        addFVec( name, value);

      // Check for occurence of an int vector and add to mvec map.
      } else if (tag == "<mvec" || tag == "<mvecfix") {
        vector<int> value = intVectorAttributeValue( line, "default=");
        int minVal = intAttributeValue( line, "min=");
        int maxVal = intAttributeValue( line, "max=");
        addMVec( name, value, hasMin, hasMax, minVal, maxVal);

      // Check for occurence of a double vector and add to pvec map.
      } else if (tag == "<pvec" || tag == "<pvecfix") {
        vector<double> value = doubleVectorAttributeValue( line, "default=");
        double minVal = doubleAttributeValue( line, "min=");
        double maxVal = doubleAttributeValue( line, "max=");
        addPVec( name, value, hasMin, hasMax, minVal, maxVal);
      }

    // End of loop over lines in input file and loop over files.
    };
  };

  // Set up default e+e- and pp tunes, if positive.
  int eeTune = mode("Tune:ee");
  if (eeTune > 0) initTuneEE( eeTune);
  int ppTune = mode("Tune:pp");
  if (ppTune > 0) initTunePP( ppTune);

  // Done.
  if (nError > 0) return false;
  isInit = true;
  return true;

}

//--------------------------------------------------------------------------

// Overwrite existing database by reading from specific file.

bool Settings::reInit(string startFile, ostream& os) {

  // Reset maps to empty.
  flags.clear();
  modes.clear();
  parms.clear();
  words.clear();
  fvecs.clear();
  mvecs.clear();
  pvecs.clear();

  // Then let normal init do the rest.
  isInit = false;
  return init(startFile, false, os);

}

//--------------------------------------------------------------------------

// Read in updates from a character string, like a line of a file.
// Is used by readString (and readFile) in Pythia.

bool Settings::readString(string line, bool warn, ostream& os) {

  // If empty line then done.
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return true;

  // If first character is not a letter, then taken to be a comment line.
  string lineNow = line;
  int firstChar = lineNow.find_first_not_of(" \n\t\v\b\r\f\a");
  if (!isalpha(lineNow[firstChar])) return true;

  // Replace an equal sign by a blank to make parsing simpler.
  while (lineNow.find("=") != string::npos) {
    int firstEqual = lineNow.find_first_of("=");
    lineNow.replace(firstEqual, 1, " ");
  }

  // Get first word of a line.
  istringstream splitLine(lineNow);
  string name;
  splitLine >> name;

  // Replace two colons by one (:: -> :) to allow for such mistakes.
  while (name.find("::") != string::npos) {
    int firstColonColon = name.find_first_of("::");
    name.replace(firstColonColon, 2, ":");
  }

  // Check whether this is in the database.
  int inDataBase = 0;
  if      (isFlag(name)) inDataBase = 1;
  else if (isMode(name)) inDataBase = 2;
  else if (isParm(name)) inDataBase = 3;
  else if (isWord(name)) inDataBase = 4;
  else if (isFVec(name)) inDataBase = 5;
  else if (isMVec(name)) inDataBase = 6;
  else if (isPVec(name)) inDataBase = 7;

  // For backwards compatibility: old (parts of) names mapped onto new ones.
  // This code currently has no use, but is partly preserved for the day
  // it may be needed again.
  /*
  if (inDataBase == 0) {
    bool retry = false;
    string nameLower = toLower(name);
    if (!retry && nameLower.find("minbias") != string::npos) {
      int firstMB = nameLower.find_first_of("minbias");
      name.replace(firstMB, 7, "nonDiffractive");
      retry = true;
    }
    if (retry) {
      if      (isFlag(name)) inDataBase = 1;
      else if (isMode(name)) inDataBase = 2;
      else if (isParm(name)) inDataBase = 3;
      else if (isWord(name)) inDataBase = 4;
      else if (isFVec(name)) inDataBase = 5;
      else if (isMVec(name)) inDataBase = 6;
      else if (isPVec(name)) inDataBase = 7;
    }
  }
  */

  // Warn and done if not in database.
  if (inDataBase == 0) {
    if (warn) os << "\n PYTHIA Error: input string not found in settings"
      << " databases::\n   " << line << endl;
    readingFailedSave = true;
    return false;
  }

  // Find value. Warn if none found.
  string valueString;
  splitLine >> valueString;
  if (!splitLine) {
    if (warn) os << "\n PYTHIA Error: variable recognized, but its value"
      << " not meaningful:\n   " << line << endl;
    readingFailedSave = true;
    return false;
  }

  // Update flag map; allow many ways to say yes.
  if (inDataBase == 1) {
    bool value = boolString(valueString);
    flag(name, value);

  // Update mode map.
  } else if (inDataBase == 2) {
    istringstream modeData(valueString);
    int value;
    modeData >> value;
    if (!modeData) {
      if (warn) os << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    if (!mode(name, value)) {
      if (warn) os << "\n PYTHIA Error: variable recognized, but its value"
        << " non-existing option:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }

  // Update parm map.
  } else if (inDataBase == 3) {
    istringstream parmData(valueString);
    double value;
    parmData >> value;
    if (!parmData) {
      if (warn) os << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    parm(name, value);

  // Update word map.
  } else if (inDataBase == 4)  {
    word(name, valueString);

  // Update fvec map.
  } else if (inDataBase == 5) {
    istringstream fvecData(valueString);
    vector<bool> value(boolVectorAttributeValue(
      "value=\"" + valueString + "\"", "value="));
    if (!fvecData) {
      if (warn) os << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    fvec(name, value);

  // Update mvec map.
  } else if (inDataBase == 6) {
    istringstream mvecData(valueString);
    vector<int> value(intVectorAttributeValue(
      "value=\"" + valueString + "\"", "value="));
    if (!mvecData) {
      if (warn) os << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    mvec(name, value);

  // Update pvec map.
  } else if (inDataBase == 7) {
    istringstream pvecData(valueString);
    vector<double> value(doubleVectorAttributeValue(
      "value=\"" + valueString + "\"", "value="));
    if (!pvecData) {
      if (warn) os << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    pvec(name, value);
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Write updates or everything to user-defined file.

bool Settings::writeFile(string toFile, bool writeAll) {

  // Open file for writing.
  const char* cstring = toFile.c_str();
  ofstream os(cstring);
  if (!os) {
    infoPtr->errorMsg("Error in Settings::writeFile:"
      " could not open file", toFile);
    return false;
  }

  // Hand over real work to next method.
  return writeFile( os, writeAll);

}

//--------------------------------------------------------------------------

// Write updates or everything to user-defined stream (or file).

bool Settings::writeFile(ostream& os, bool writeAll) {

  // Write simple header as comment.
  if (writeAll) os << "! List of all current PYTHIA ";
  else          os << "! List of all modified PYTHIA ";
  os << fixed << setprecision(3) << parm("Pythia:versionNumber")
     << " settings.\n";

  // Iterators for the flag, mode and parm tables.
  map<string, Flag>::iterator flagEntry = flags.begin();
  map<string, Mode>::iterator modeEntry = modes.begin();
  map<string, Parm>::iterator parmEntry = parms.begin();
  map<string, Word>::iterator wordEntry = words.begin();
  map<string, FVec>::iterator fvecEntry = fvecs.begin();
  map<string, MVec>::iterator mvecEntry = mvecs.begin();
  map<string, PVec>::iterator pvecEntry = pvecs.begin();

  // Loop while there is something left to do.
  while (flagEntry != flags.end() || modeEntry != modes.end()
      || parmEntry != parms.end() || wordEntry != words.end()
      || fvecEntry != fvecs.end() || mvecEntry != mvecs.end()
      || pvecEntry != pvecs.end() ) {

    // Check if a flag is next in lexigraphical order; if so print it.
    if ( flagEntry != flags.end()
      && ( modeEntry == modes.end() || flagEntry->first < modeEntry->first )
      && ( parmEntry == parms.end() || flagEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || flagEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || flagEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || flagEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || flagEntry->first < pvecEntry->first )
      ) {
      string state[2] = {"off", "on"};
      bool valNow = flagEntry->second.valNow;
      bool valDefault = flagEntry->second.valDefault;
      if ( writeAll || valNow != valDefault )
        os << flagEntry->second.name << " = " << state[valNow] << "\n";
      ++flagEntry;

    // Else check if mode is next, and if so print it.
    } else if ( modeEntry != modes.end()
      && ( parmEntry == parms.end() || modeEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || modeEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || modeEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || modeEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || modeEntry->first < pvecEntry->first )
      ) {
      int valNow = modeEntry->second.valNow;
      int valDefault = modeEntry->second.valDefault;
      if ( writeAll || valNow != valDefault )
        os << modeEntry->second.name << " = " << valNow << "\n";
      ++modeEntry;

    // Else check if parm is next, and if so print it;
    // fixed or scientific depending on value.
    } else if ( parmEntry != parms.end()
      && ( wordEntry == words.end() || parmEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || parmEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || parmEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || parmEntry->first < pvecEntry->first )
      ) {
      double valNow = parmEntry->second.valNow;
      double valDefault = parmEntry->second.valDefault;
      if ( writeAll || valNow != valDefault ) {
        os  << parmEntry->second.name << " = ";
        if ( valNow == 0. ) os << fixed << setprecision(1);
        else if ( abs(valNow) < 0.001 ) os << scientific << setprecision(4);
        else if ( abs(valNow) < 0.1 ) os << fixed << setprecision(7);
        else if ( abs(valNow) < 1000. ) os << fixed << setprecision(5);
        else if ( abs(valNow) < 1000000. ) os << fixed << setprecision(3);
        else os << scientific << setprecision(4);
        os << valNow << "\n";
      }
      ++parmEntry;

    // Else check if word is next, and if so print it.
    } else  if ( wordEntry != words.end()
      && ( fvecEntry == fvecs.end() || wordEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || wordEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || wordEntry->first < pvecEntry->first )
      ) {
      string valNow = wordEntry->second.valNow;
      string valDefault = wordEntry->second.valDefault;
      if ( writeAll || valNow != valDefault )
        os << wordEntry->second.name << " = " << valNow << "\n";
      ++wordEntry;

    // Else check if fvec is next, and if so print it.
    } else if ( fvecEntry != fvecs.end()
      && ( mvecEntry == mvecs.end() || fvecEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || fvecEntry->first < pvecEntry->first )
      ) {
      string state[2] = {"off", "on"};
      vector<bool> valNow = fvecEntry->second.valNow;
      vector<bool> valDefault = fvecEntry->second.valDefault;
      if ( writeAll || valNow != valDefault ) {
        os  << fvecEntry->second.name << " = ";
        for (vector<bool>::iterator val = valNow.begin();
             val != --valNow.end(); ++val) os << state[*val] << ",";
        os << *(--valNow.end()) << "\n";
      }
      ++fvecEntry;

    // Else check if mvec is next, and if so print it.
    } else if ( mvecEntry != mvecs.end()
      && ( pvecEntry == pvecs.end() || mvecEntry->first < pvecEntry->first )
      ) {
      vector<int> valNow = mvecEntry->second.valNow;
      vector<int> valDefault = mvecEntry->second.valDefault;
      if ( writeAll || valNow != valDefault ) {
        os  << mvecEntry->second.name << " = ";
        for (vector<int>::iterator val = valNow.begin();
             val != --valNow.end(); ++val) os << *val << ",";
        os << *(--valNow.end()) << "\n";
      }
      ++mvecEntry;

    // Else print pvec; fixed or scientific depending on value.
    } else {
      vector<double> valNow = pvecEntry->second.valNow;
      vector<double> valDefault = pvecEntry->second.valDefault;
      if ( writeAll || valNow != valDefault ) {
        os  << pvecEntry->second.name << " = ";
        for (vector<double>::iterator val = valNow.begin();
             val != --valNow.end(); ++val) {
          if ( *val == 0. ) os << fixed << setprecision(1);
          else if ( abs(*val) < 0.001 ) os << scientific << setprecision(4);
          else if ( abs(*val) < 0.1 ) os << fixed << setprecision(7);
          else if ( abs(*val) < 1000. ) os << fixed << setprecision(5);
          else if ( abs(*val) < 1000000. ) os << fixed << setprecision(3);
          else os << scientific << setprecision(4);
          os << *val << ",";
        } os << *(--valNow.end()) << "\n";
      }
      ++pvecEntry;
    }
  } ;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Print out table of database in lexigraphical order.

void Settings::list(bool doListAll,  bool doListString, string match,
  ostream& os) {

  // Table header; output for bool as off/on.
  if (doListAll)
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Word + FVec + MVec + PVec "
       << "Settings (all)  ----------------------------------* \n";
  else if (!doListString)
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Word + FVec + MVec + PVec "
       << "Settings (changes only)  -------------------------* \n" ;
  else
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Word + FVec + MVec + PVec "
       << "Settings (with requested string) -----------------* \n" ;
  os << " |                                                           "
     << "                                                      | \n"
     << " | Name                                          |           "
     << "           Now |      Default         Min         Max | \n"
     << " |                                               |           "
     << "               |                                      | \n";

  // Convert input string to lowercase for match.
  match = toLower(match);
  if (match == "") match = "             ";

  // Iterators for the flag, mode and parm tables.
  map<string, Flag>::iterator flagEntry = flags.begin();
  map<string, Mode>::iterator modeEntry = modes.begin();
  map<string, Parm>::iterator parmEntry = parms.begin();
  map<string, Word>::iterator wordEntry = words.begin();
  map<string, FVec>::iterator fvecEntry = fvecs.begin();
  map<string, MVec>::iterator mvecEntry = mvecs.begin();
  map<string, PVec>::iterator pvecEntry = pvecs.begin();

  // Loop while there is something left to do.
  while (flagEntry != flags.end() || modeEntry != modes.end()
      || parmEntry != parms.end() || wordEntry != words.end()
      || fvecEntry != fvecs.end() || mvecEntry != mvecs.end()
      || pvecEntry != pvecs.end() ) {

    // Check if a flag is next in lexigraphical order; if so print it.
    if ( flagEntry != flags.end()
      && ( modeEntry == modes.end() || flagEntry->first < modeEntry->first )
      && ( parmEntry == parms.end() || flagEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || flagEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || flagEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || flagEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || flagEntry->first < pvecEntry->first )
      ) {
      string state[2] = {"off", "on"};
      bool valNow = flagEntry->second.valNow;
      bool valDefault = flagEntry->second.valDefault;
      if ( doListAll || (!doListString && valNow != valDefault)
        || (doListString && flagEntry->first.find(match) != string::npos) )
        os << " | " << setw(45) << left
           << flagEntry->second.name << " | " << setw(24) << right
           << state[valNow] << " | " << setw(12) << state[valDefault]
           << "                         | \n";
      ++flagEntry;

    // Else check if mode is next, and if so print it.
    } else if ( modeEntry != modes.end()
      && ( parmEntry == parms.end() || modeEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || modeEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || modeEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || modeEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || modeEntry->first < pvecEntry->first )
      ) {
      int valNow = modeEntry->second.valNow;
      int valDefault = modeEntry->second.valDefault;
      if ( doListAll || (!doListString && valNow != valDefault)
        || (doListString && modeEntry->first.find(match) != string::npos) ) {
        os << " | " << setw(45) << left
           << modeEntry->second.name << " | " << setw(24) << right
           << valNow << " | " << setw(12) << valDefault;
        if (modeEntry->second.hasMin)
          os << setw(12) << modeEntry->second.valMin;
        else os << "            ";
        if (modeEntry->second.hasMax)
          os << setw(12) << modeEntry->second.valMax;
        else os << "            ";
        os << " | \n";
      }
      ++modeEntry;

    // Else check if parm is next, and if so print it;
    // fixed or scientific depending on value.
    } else if ( parmEntry != parms.end()
      && ( wordEntry == words.end() || parmEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || parmEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || parmEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || parmEntry->first < pvecEntry->first )
      ) {
      double valNow = parmEntry->second.valNow;
      double valDefault = parmEntry->second.valDefault;
      if ( doListAll || (!doListString && valNow != valDefault )
        || (doListString && parmEntry->first.find(match) != string::npos) ) {
        os << " | " << setw(45) << left
           << parmEntry->second.name << right << " |             ";
        for (int i = 0; i < 4; ++i) {
          if (i == 1) valNow = valDefault;
          if (i == 2) valNow = parmEntry->second.valMin;
          if (i == 3) valNow = parmEntry->second.valMax;
          if ( (i == 2 && !parmEntry->second.hasMin)
            || (i == 3 && !parmEntry->second.hasMax) )
            os << "            ";
          else if ( valNow == 0. )
            os << fixed << setprecision(1) << setw(12) << valNow;
          else if ( abs(valNow) < 0.001 )
            os << scientific << setprecision(4) << setw(12) << valNow;
          else if ( abs(valNow) < 0.1 )
            os << fixed << setprecision(7) << setw(12) << valNow;
          else if ( abs(valNow) < 1000. )
            os << fixed << setprecision(5) << setw(12) << valNow;
          else if ( abs(valNow) < 1000000. )
            os << fixed << setprecision(3) << setw(12) << valNow;
          else
            os << scientific << setprecision(4) << setw(12) << valNow;
          if (i == 0) os << " | ";
        }
        os << " | \n";
      }
      ++parmEntry;

    // Else check if word is next, and if so print it.
    } else  if ( wordEntry != words.end()
      && ( fvecEntry == fvecs.end() || wordEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || wordEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || wordEntry->first < pvecEntry->first )
      ) {
      string valNow = wordEntry->second.valNow;
      string valDefault = wordEntry->second.valDefault;
      int blankLeft = max(0, 60 - max(24, int(valNow.length()) )
        - max(12, int(valDefault.length()) ) );
      string blankPad( blankLeft, ' ');
      if ( doListAll || (!doListString && valNow != valDefault)
        || (doListString && wordEntry->first.find(match) != string::npos) )
        os << " | " << setw(45) << left
           << wordEntry->second.name << " | " << setw(24) << right
           << valNow << " | " << setw(12) << valDefault << blankPad
           << " | \n";
      ++wordEntry;

    // Else check if fvec is next, and if so print it.
    } else if ( fvecEntry != fvecs.end()
      && ( mvecEntry == mvecs.end() || fvecEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || fvecEntry->first < pvecEntry->first )
      ) {
      string state[2] = {"off", "on"};
      vector<bool> valsNow = fvecEntry->second.valNow;
      vector<bool> valsDefault = fvecEntry->second.valDefault;
      bool valNow(false), valDefault(false);
      if ( doListAll || (!doListString && valsNow != valsDefault )
        || (doListString && fvecEntry->first.find(match) != string::npos) ) {
        for (unsigned int i = 0; i < valsNow.size() || i < valsDefault.size();
             ++i) {
          if ( i == 0 )
            os << " | " << setw(45) << left
               << fvecEntry->second.name << right << " |             ";
          else
            os << " | " << setw(45) << " " << right << " |             ";
          for (int j = 0; j < 4; ++j) {
            if (i < valsNow.size()) valNow = valsNow[i];
            if (i < valsDefault.size()) valDefault = valsDefault[i];
            if (j == 1) valNow = valDefault;
            if ( (j == 0 && i >= valsNow.size())
                 || (j == 1 && i >= valsDefault.size()) || (j > 1) )
              os << "            ";
            else os << setw(12) << state[valNow];
            if (j == 0) os << " | ";
          }
          os << " | \n";
        }
      }
      ++fvecEntry;

    // Else check if mvec is next, and if so print it.
    } else if ( mvecEntry != mvecs.end()
      && ( pvecEntry == pvecs.end() || mvecEntry->first < pvecEntry->first )
      ) {
      vector<int> valsNow = mvecEntry->second.valNow;
      vector<int> valsDefault = mvecEntry->second.valDefault;
      int valNow(0), valDefault(0);
      if ( doListAll || (!doListString && valsNow != valsDefault )
        || (doListString && mvecEntry->first.find(match) != string::npos) ) {
        for (unsigned int i = 0; i < valsNow.size() || i < valsDefault.size();
             ++i) {
          if ( i == 0 )
            os << " | " << setw(45) << left
               << mvecEntry->second.name << right << " |             ";
          else
            os << " | " << setw(45) << " " << right << " |             ";
          for (int j = 0; j < 4; ++j) {
            if (i < valsNow.size()) valNow = valsNow[i];
            if (i < valsDefault.size()) valDefault = valsDefault[i];
            if (j == 1) valNow = valDefault;
            if (j == 2) valNow = mvecEntry->second.valMin;
            if (j == 3) valNow = mvecEntry->second.valMax;
            if ( (j == 0 && i >= valsNow.size())
                 || (j == 1 && i >= valsDefault.size())
                 || (j == 2 && !mvecEntry->second.hasMin)
                 || (j == 3 && !mvecEntry->second.hasMax) )
              os << "            ";
            else os << setw(12) << valNow;
            if (j == 0) os << " | ";
          }
          os << " | \n";
        }
      }
      ++mvecEntry;

    // Else print pvec; fixed or scientific depending on value.
    } else {
      vector<double> valsNow = pvecEntry->second.valNow;
      vector<double> valsDefault = pvecEntry->second.valDefault;
      double valNow(0), valDefault(0);
      if ( doListAll || (!doListString && valsNow != valsDefault )
        || (doListString && pvecEntry->first.find(match) != string::npos) ) {
        for (unsigned int i = 0; i < valsNow.size() || i < valsDefault.size();
             ++i) {
          if ( i == 0 )
            os << " | " << setw(45) << left
               << pvecEntry->second.name << right << " |             ";
          else
            os << " | " << setw(45) << " " << right << " |             ";
          for (int j = 0; j < 4; ++j) {
            if (i < valsNow.size()) valNow = valsNow[i];
            if (i < valsDefault.size()) valDefault = valsDefault[i];
            if (j == 1) valNow = valDefault;
            if (j == 2) valNow = pvecEntry->second.valMin;
            if (j == 3) valNow = pvecEntry->second.valMax;
            if ( (j == 0 && i >= valsNow.size())
                 || (j == 1 && i >= valsDefault.size())
                 || (j == 2 && !pvecEntry->second.hasMin)
                 || (j == 3 && !pvecEntry->second.hasMax) )
              os << "            ";
            else if ( valNow == 0. )
              os << fixed << setprecision(1) << setw(12) << valNow;
            else if ( abs(valNow) < 0.001 )
              os << scientific << setprecision(4) << setw(12) << valNow;
            else if ( abs(valNow) < 0.1 )
              os << fixed << setprecision(7) << setw(12) << valNow;
            else if ( abs(valNow) < 1000. )
              os << fixed << setprecision(5) << setw(12) << valNow;
            else if ( abs(valNow) < 1000000. )
              os << fixed << setprecision(3) << setw(12) << valNow;
            else
              os << scientific << setprecision(4) << setw(12) << valNow;
            if (j == 0) os << " | ";
          }
          os << " | \n";
        }
      }
      ++pvecEntry;

    }
  } ;

  // End of loop over database contents.
  os << " |                                                           "
     << "                                                      | \n"
     << " *-------  End PYTHIA Flag + Mode + Parm + Word + FVec + MVec + PVec "
     << "Settings  ------------------------------------* " << endl;

}

//--------------------------------------------------------------------------

// Reset all values to their defaults.

void Settings::resetAll() {

  // Loop through the flags table, resetting all entries.
  for (map<string, Flag>::iterator flagEntry = flags.begin();
    flagEntry != flags.end(); ++flagEntry) {
    string name = flagEntry->first;
    resetFlag(name);
  }

  // Loop through the modes table, resetting all entries.
  for (map<string, Mode>::iterator modeEntry = modes.begin();
    modeEntry != modes.end(); ++modeEntry) {
    string name = modeEntry->first;
    resetMode(name);
  }

  // Loop through the parms table, resetting all entries.
  for (map<string, Parm>::iterator parmEntry = parms.begin();
    parmEntry != parms.end(); ++parmEntry) {
    string name = parmEntry->first;
    resetParm(name);
  }

  // Loop through the words table, resetting all entries.
  for (map<string, Word>::iterator wordEntry = words.begin();
    wordEntry != words.end(); ++wordEntry) {
    string name = wordEntry->first;
    resetWord(name);
  }

  // Loop through the fvecs table, resetting all entries.
  for (map<string, FVec>::iterator fvecEntry = fvecs.begin();
    fvecEntry != fvecs.end(); ++fvecEntry) {
    string name = fvecEntry->first;
    resetFVec(name);
  }

  // Loop through the mvecs table, resetting all entries.
  for (map<string, MVec>::iterator mvecEntry = mvecs.begin();
    mvecEntry != mvecs.end(); ++mvecEntry) {
    string name = mvecEntry->first;
    resetMVec(name);
  }

  // Loop through the pvecs table, resetting all entries.
  for (map<string, PVec>::iterator pvecEntry = pvecs.begin();
    pvecEntry != pvecs.end(); ++pvecEntry) {
    string name = pvecEntry->first;
    resetPVec(name);
  }

}

//--------------------------------------------------------------------------

// Give back current value, with check that key exists.

bool Settings::flag(string keyIn) {
  if (isFlag(keyIn)) return flags[toLower(keyIn)].valNow;
  infoPtr->errorMsg("Error in Settings::flag: unknown key", keyIn);
  return false;
}

int Settings::mode(string keyIn) {
  if (isMode(keyIn)) return modes[toLower(keyIn)].valNow;
  infoPtr->errorMsg("Error in Settings::mode: unknown key", keyIn);
  return 0;
}

double Settings::parm(string keyIn) {
  if (isParm(keyIn)) return parms[toLower(keyIn)].valNow;
  infoPtr->errorMsg("Error in Settings::parm: unknown key", keyIn);
  return 0.;
}

string Settings::word(string keyIn) {
  if (isWord(keyIn)) return words[toLower(keyIn)].valNow;
  infoPtr->errorMsg("Error in Settings::word: unknown key", keyIn);
  return " ";
}

vector<bool> Settings::fvec(string keyIn) {
  if (isFVec(keyIn)) return fvecs[toLower(keyIn)].valNow;
  infoPtr->errorMsg("Error in Settings::fvec: unknown key", keyIn);
  return vector<bool>(1, false);
}

vector<int> Settings::mvec(string keyIn) {
  if (isMVec(keyIn)) return mvecs[toLower(keyIn)].valNow;
  infoPtr->errorMsg("Error in Settings::mvec: unknown key", keyIn);
  return vector<int>(1, 0);
}

vector<double> Settings::pvec(string keyIn) {
  if (isPVec(keyIn)) return pvecs[toLower(keyIn)].valNow;
  infoPtr->errorMsg("Error in Settings::pvec: unknown key", keyIn);
  return vector<double>(1, 0.);
}

//--------------------------------------------------------------------------

// Give back default value, with check that key exists.

bool Settings::flagDefault(string keyIn) {
  if (isFlag(keyIn)) return flags[toLower(keyIn)].valDefault;
  infoPtr->errorMsg("Error in Settings::flagDefault: unknown key", keyIn);
  return false;
}

int Settings::modeDefault(string keyIn) {
  if (isMode(keyIn)) return modes[toLower(keyIn)].valDefault;
  infoPtr->errorMsg("Error in Settings::modeDefault: unknown key", keyIn);
  return 0;
}

double Settings::parmDefault(string keyIn) {
  if (isParm(keyIn)) return parms[toLower(keyIn)].valDefault;
  infoPtr->errorMsg("Error in Settings::parmDefault: unknown key", keyIn);
  return 0.;
}

string Settings::wordDefault(string keyIn) {
  if (isWord(keyIn)) return words[toLower(keyIn)].valDefault;
  infoPtr->errorMsg("Error in Settings::wordDefault: unknown key", keyIn);
  return " ";
}

vector<bool> Settings::fvecDefault(string keyIn) {
  if (isFVec(keyIn)) return fvecs[toLower(keyIn)].valDefault;
  infoPtr->errorMsg("Error in Settings::fvecDefault: unknown key", keyIn);
  return vector<bool>(1, false);
}

vector<int> Settings::mvecDefault(string keyIn) {
  if (isMVec(keyIn)) return mvecs[toLower(keyIn)].valDefault;
  infoPtr->errorMsg("Error in Settings::mvecDefault: unknown key", keyIn);
  return vector<int>(1, 0);
}

vector<double> Settings::pvecDefault(string keyIn) {
  if (isPVec(keyIn)) return pvecs[toLower(keyIn)].valDefault;
  infoPtr->errorMsg("Error in Settings::pvecDefault: unknown key", keyIn);
  return vector<double>(1, 0.);
}

//--------------------------------------------------------------------------

// Get a map of entries whose names contain the string "match".

map<string, Flag> Settings::getFlagMap(string match) {
  // Make the match string lower case. Start with an empty map.
  match = toLower(match);
  map<string, Flag> flagMap;
  // Loop over the flag map (using iterator).
  for (map<string,Flag>::iterator flagEntry = flags.begin();
       flagEntry != flags.end(); ++flagEntry)
    if (flagEntry->first.find(match) != string::npos)
      flagMap[flagEntry->first] = flagEntry->second;
  return flagMap;
}

map<string, Mode> Settings::getModeMap(string match) {
  // Make the match string lower case. Start with an empty map.
  match = toLower(match);
  map<string, Mode> modeMap;
  // Loop over the mode map (using iterator).
  for (map<string,Mode>::iterator modeEntry = modes.begin();
       modeEntry != modes.end(); ++modeEntry)
    if (modeEntry->first.find(match) != string::npos)
      modeMap[modeEntry->first] = modeEntry->second;
  return modeMap;
}

map<string, Parm> Settings::getParmMap(string match) {
  // Make the match string lower case. Start with an empty map.
  match = toLower(match);
  map<string, Parm> parmMap;
  // Loop over the parm map (using iterator).
  for (map<string,Parm>::iterator parmEntry = parms.begin();
       parmEntry != parms.end(); ++parmEntry)
    if (parmEntry->first.find(match) != string::npos)
      parmMap[parmEntry->first] = parmEntry->second;
  return parmMap;
}

map<string, Word> Settings::getWordMap(string match) {
  // Make the match string lower case. Start with an empty map.
  match = toLower(match);
  map<string, Word> wordMap;
  // Loop over the word map (using iterator).
  for (map<string,Word>::iterator wordEntry = words.begin();
       wordEntry != words.end(); ++wordEntry)
    if (wordEntry->first.find(match) != string::npos)
      wordMap[wordEntry->first] = wordEntry->second;
  return wordMap;
}

map<string, FVec> Settings::getFVecMap(string match) {
  // Make the match string lower case. Start with an empty map.
  match = toLower(match);
  map<string, FVec> fvecMap;
  // Loop over the fvec map (using iterator).
  for (map<string,FVec>::iterator fvecEntry = fvecs.begin();
       fvecEntry != fvecs.end(); ++fvecEntry)
    if (fvecEntry->first.find(match) != string::npos)
      fvecMap[fvecEntry->first] = fvecEntry->second;
  return fvecMap;
}

map<string, MVec> Settings::getMVecMap(string match) {
  // Make the match string lower case. Start with an empty map.
  match = toLower(match);
  map<string, MVec> mvecMap;
  // Loop over the mvec map (using iterator).
  for (map<string,MVec>::iterator mvecEntry = mvecs.begin();
       mvecEntry != mvecs.end(); ++mvecEntry)
    if (mvecEntry->first.find(match) != string::npos)
      mvecMap[mvecEntry->first] = mvecEntry->second;
  return mvecMap;
}

map<string, PVec> Settings::getPVecMap(string match) {
  // Make the match string lower case. Start with an empty map.
  match = toLower(match);
  map<string, PVec> pvecMap;
  // Loop over the pvec map (using iterator).
  for (map<string,PVec>::iterator pvecEntry = pvecs.begin();
       pvecEntry != pvecs.end(); ++pvecEntry)
    if (pvecEntry->first.find(match) != string::npos)
      pvecMap[pvecEntry->first] = pvecEntry->second;
  return pvecMap;
}

//--------------------------------------------------------------------------

// Change current value, respecting limits.

void Settings::flag(string keyIn, bool nowIn) {
  string keyLower = toLower(keyIn);
  if (isFlag(keyIn)) flags[keyLower].valNow = nowIn;
  // Print:quiet  triggers a whole set of changes.
  if (keyLower == "print:quiet") printQuiet( nowIn);
}

bool Settings:: mode(string keyIn, int nowIn) {
  if (isMode(keyIn)) {
    string keyLower = toLower(keyIn);
    Mode& modeNow = modes[keyLower];
    // For modepick and modefix fail if values are outside range.
    if (modeNow.optOnly && (nowIn < modeNow.valMin || nowIn > modeNow.valMax) )
      return false;
    if (modeNow.hasMin && nowIn < modeNow.valMin)
      modeNow.valNow = modeNow.valMin;
    else if (modeNow.hasMax && nowIn > modeNow.valMax)
      modeNow.valNow = modeNow.valMax;
    else modeNow.valNow = nowIn;
    // Tune:ee and Tune:pp each trigger a whole set of changes.
    if (keyLower == "tune:ee") initTuneEE( modeNow.valNow);
    if (keyLower == "tune:pp") initTunePP( modeNow.valNow);
  }
  return true;

}

void Settings::parm(string keyIn, double nowIn) {
  if (isParm(keyIn)) {
    Parm& parmNow = parms[toLower(keyIn)];
    if (parmNow.hasMin && nowIn < parmNow.valMin)
      parmNow.valNow = parmNow.valMin;
    else if (parmNow.hasMax && nowIn > parmNow.valMax)
      parmNow.valNow = parmNow.valMax;
    else parmNow.valNow = nowIn;
  }
}

void Settings::word(string keyIn, string nowIn) {
    if (isWord(keyIn)) words[toLower(keyIn)].valNow = nowIn;
}

void Settings::fvec(string keyIn, vector<bool> nowIn) {
  if (isFVec(keyIn)) {
    FVec& fvecNow = fvecs[toLower(keyIn)];
    fvecNow.valNow.clear();
    for(vector<bool>::iterator now = nowIn.begin();
        now != nowIn.end(); now++)
      fvecNow.valNow.push_back(*now);
  }
}

void Settings::mvec(string keyIn, vector<int> nowIn) {
  if (isMVec(keyIn)) {
    MVec& mvecNow = mvecs[toLower(keyIn)];
    mvecNow.valNow.clear();
    for(vector<int>::iterator now = nowIn.begin();
        now != nowIn.end(); now++) {
      if (mvecNow.hasMin && *now < mvecNow.valMin)
        mvecNow.valNow.push_back(mvecNow.valMin);
      else if (mvecNow.hasMax && *now > mvecNow.valMax)
        mvecNow.valNow.push_back(mvecNow.valMax);
      else mvecNow.valNow.push_back(*now);
    }
  }
}

void Settings::pvec(string keyIn, vector<double> nowIn) {
  if (isPVec(keyIn)) {
    PVec& pvecNow = pvecs[toLower(keyIn)];
    pvecNow.valNow.clear();
    for(vector<double>::iterator now = nowIn.begin();
        now != nowIn.end(); now++) {
      if (pvecNow.hasMin && *now < pvecNow.valMin)
        pvecNow.valNow.push_back(pvecNow.valMin);
      else if (pvecNow.hasMax && *now > pvecNow.valMax)
        pvecNow.valNow.push_back(pvecNow.valMax);
      else pvecNow.valNow.push_back(*now);
    }
  }
}

//--------------------------------------------------------------------------

// Change current value, disregarding limits.

void Settings::forceMode(string keyIn, int nowIn) {
  if (isMode(keyIn)) {
    string keyLower = toLower(keyIn);
    Mode& modeNow   = modes[keyLower];
    modeNow.valNow  = nowIn;
    // Tune:ee and Tune:pp each trigger a whole set of changes.
    if (keyLower == "tune:ee") initTuneEE( modeNow.valNow);
    if (keyLower == "tune:pp") initTunePP( modeNow.valNow);
  }
}

void Settings::forceParm(string keyIn, double nowIn) {
  if (isParm(keyIn)) parms[toLower(keyIn)].valNow = nowIn;
}

void Settings::forceMVec(string keyIn, vector<int> nowIn) {
  if (isMVec(keyIn)) mvecs[toLower(keyIn)].valNow = nowIn;
}

void Settings::forcePVec(string keyIn, vector<double> nowIn) {
  if (isPVec(keyIn)) pvecs[toLower(keyIn)].valNow = nowIn;
}

//--------------------------------------------------------------------------

// Restore current value to default.

void Settings::resetFlag(string keyIn) {
  if (isFlag(keyIn)) flags[toLower(keyIn)].valNow
    = flags[toLower(keyIn)].valDefault ;
}

void Settings::resetMode(string keyIn) {
  string keyLower = toLower(keyIn);
  if (isMode(keyIn)) modes[keyLower].valNow
    = modes[toLower(keyIn)].valDefault ;
  // For Tune:ee and Tune:pp must also restore variables involved in tunes.
  if (keyLower == "tune:ee") resetTuneEE();
  if (keyLower == "tune:pp") resetTunePP();
}

void Settings::resetParm(string keyIn) {
  if (isParm(keyIn)) parms[toLower(keyIn)].valNow
    = parms[toLower(keyIn)].valDefault ;
}

void Settings::resetWord(string keyIn) {
  if (isWord(keyIn)) words[toLower(keyIn)].valNow
    = words[toLower(keyIn)].valDefault ;
}

void Settings::resetFVec(string keyIn) {
  if (isFVec(keyIn)) fvecs[toLower(keyIn)].valNow
    = fvecs[toLower(keyIn)].valDefault ;
}

void Settings::resetMVec(string keyIn) {
  if (isMVec(keyIn)) mvecs[toLower(keyIn)].valNow
    = mvecs[toLower(keyIn)].valDefault ;
}

void Settings::resetPVec(string keyIn) {
  if (isPVec(keyIn)) pvecs[toLower(keyIn)].valNow
    = pvecs[toLower(keyIn)].valDefault ;
}

//--------------------------------------------------------------------------

// Regulate level of printout by overall change of settings.

void Settings::printQuiet(bool quiet) {

  // Switch off as much output as possible.
  if (quiet) {
    flag("Init:showProcesses",               false );
    flag("Init:showMultipartonInteractions", false );
    flag("Init:showChangedSettings",         false );
    flag("Init:showAllSettings",             false );
    flag("Init:showChangedParticleData",     false );
    flag("Init:showChangedResonanceData",    false );
    flag("Init:showAllParticleData",         false );
    mode("Init:showOneParticleData",             0 );
    mode("Next:numberCount",                     0 );
    mode("Next:numberShowLHA",                   0 );
    mode("Next:numberShowInfo",                  0 );
    mode("Next:numberShowProcess",               0 );
    mode("Next:numberShowEvent",                 0 );

  // Restore ouput settings to default.
  } else {
    resetFlag("Init:showProcesses");
    resetFlag("Init:showMultipartonInteractions");
    resetFlag("Init:showChangedSettings");
    resetFlag("Init:showAllSettings");
    resetFlag("Init:showChangedParticleData");
    resetFlag("Init:showChangedResonanceData");
    resetFlag("Init:showAllParticleData");
    resetMode("Init:showOneParticleData");
    resetMode("Next:numberCount");
    resetMode("Next:numberShowLHA");
    resetMode("Next:numberShowInfo");
    resetMode("Next:numberShowProcess");
    resetMode("Next:numberShowEvent");
  }

}

//--------------------------------------------------------------------------

// Restore all e+e- settings to their original values.

void Settings::resetTuneEE() {

  // Flavour composition.
  resetParm("StringFlav:probStoUD");
  resetParm("StringFlav:probQQtoQ");
  resetParm("StringFlav:probSQtoQQ");
  resetParm("StringFlav:probQQ1toQQ0");
  resetParm("StringFlav:mesonUDvector");
  resetParm("StringFlav:mesonSvector");
  resetParm("StringFlav:mesonCvector");
  resetParm("StringFlav:mesonBvector");
  resetParm("StringFlav:etaSup");
  resetParm("StringFlav:etaPrimeSup");
  resetParm("StringFlav:popcornSpair");
  resetParm("StringFlav:popcornSmeson");
  resetFlag("StringFlav:suppressLeadingB");

  // String breaks: z.
  resetParm("StringZ:aLund");
  resetParm("StringZ:bLund");
  resetParm("StringZ:aExtraSquark");
  resetParm("StringZ:aExtraDiquark");
  resetParm("StringZ:rFactC");
  resetParm("StringZ:rFactB");

  // String breaks: pT.
  resetParm("StringPT:sigma");
  resetParm("StringPT:enhancedFraction");
  resetParm("StringPT:enhancedWidth");

  // FSR: strong coupling, IR cutoff.
  resetParm("TimeShower:alphaSvalue");
  resetMode("TimeShower:alphaSorder");
  resetFlag("TimeShower:alphaSuseCMW");
  resetParm("TimeShower:pTmin");
  resetParm("TimeShower:pTminChgQ");

}

//--------------------------------------------------------------------------

// Restore all pp settings to their original values.

void Settings::resetTunePP() {

  // PDF set.
  resetWord("PDF:pSet");

  // Hard matrix elements alpha_s value.
  resetParm("SigmaProcess:alphaSvalue");

  // Diffraction: cross sections and mass distributions.
  resetFlag("SigmaTotal:zeroAXB");
  resetFlag("SigmaDiffractive:dampen");
  resetParm("SigmaDiffractive:maxXB");
  resetParm("SigmaDiffractive:maxAX");
  resetParm("SigmaDiffractive:maxXX");
  resetParm("Diffraction:largeMassSuppress");

  // FSR: dipoles to beam, spin correlations.
  resetFlag("TimeShower:dampenBeamRecoil");
  resetFlag("TimeShower:phiPolAsym");

  // ISR: strong coupling, IR cutoff, coherence and spin correlations.
  resetParm("SpaceShower:alphaSvalue");
  resetMode("SpaceShower:alphaSorder");
  resetParm("SpaceShower:alphaSuseCMW");
  resetFlag("SpaceShower:samePTasMPI");
  resetParm("SpaceShower:pT0Ref");
  resetParm("SpaceShower:ecmRef");
  resetParm("SpaceShower:ecmPow");
  resetParm("SpaceShower:pTmaxFudge");
  resetParm("SpaceShower:pTdampFudge");
  resetFlag("SpaceShower:rapidityOrder");
  resetFlag("SpaceShower:phiPolAsym");
  resetFlag("SpaceShower:phiIntAsym");

  // MPI: strong coupling, IR regularization, energy scaling.
  resetParm("MultipartonInteractions:alphaSvalue");
  resetParm("MultipartonInteractions:pT0Ref");
  resetParm("MultipartonInteractions:ecmRef");
  resetParm("MultipartonInteractions:ecmPow");
  resetMode("MultipartonInteractions:bProfile");
  resetParm("MultipartonInteractions:expPow");
  resetParm("MultipartonInteractions:a1");

  // Beam remnant parameters.
  resetParm("BeamRemnants:primordialKTsoft");
  resetParm("BeamRemnants:primordialKThard");
  resetParm("BeamRemnants:halfScaleForKT");
  resetParm("BeamRemnants:halfMassForKT");

  // Colour reconnection parameters.
  resetMode("ColourReconnection:mode");
  resetParm("ColourReconnection:range");

}

//--------------------------------------------------------------------------

// Set the values related to a tune of e+e- data,
// i.e. mainly for final-state radiation and hadronization.

void Settings::initTuneEE( int eeTune) {

  // Restore all e+e- settings to their original values.
  // Is first step for setting up a specific tune.
  if (eeTune != 0) resetTuneEE();

  // Old flavour and FSR defaults carried over from very old JETSET tune,
  // only with alphaS roughly tuned for "new" pT-ordered shower.
  if (eeTune == 1) {
    parm("StringFlav:probStoUD",        0.30  );
    parm("StringFlav:probQQtoQ",        0.10  );
    parm("StringFlav:probSQtoQQ",       0.40  );
    parm("StringFlav:probQQ1toQQ0",     0.05  );
    parm("StringFlav:mesonUDvector",    1.00  );
    parm("StringFlav:mesonSvector",     1.50  );
    parm("StringFlav:mesonCvector",     2.50  );
    parm("StringFlav:mesonBvector",     3.00  );
    parm("StringFlav:etaSup",           1.00  );
    parm("StringFlav:etaPrimeSup",      0.40  );
    parm("StringFlav:popcornSpair",     0.50  );
    parm("StringFlav:popcornSmeson",    0.50  );
    flag("StringFlav:suppressLeadingB", false );
    parm("StringZ:aLund",               0.30  );
    parm("StringZ:bLund",               0.58  );
    parm("StringZ:aExtraSquark",        0.00  );
    parm("StringZ:aExtraDiquark",       0.50  );
    parm("StringZ:rFactC",              1.00  );
    parm("StringZ:rFactB",              1.00  );
    parm("StringPT:sigma",              0.36  );
    parm("StringPT:enhancedFraction",   0.01  );
    parm("StringPT:enhancedWidth",      2.0   );
    parm("TimeShower:alphaSvalue",      0.137 );
    mode("TimeShower:alphaSorder",      1     );
    flag("TimeShower:alphaSuseCMW",     false );
    parm("TimeShower:pTmin",            0.5   );
    parm("TimeShower:pTminChgQ",        0.5   );
  }

  // Marc Montull's tune to particle composition at LEP1 (August 2007).
  else if (eeTune == 2) {
    parm("StringFlav:probStoUD",        0.22  );
    parm("StringFlav:probQQtoQ",        0.08  );
    parm("StringFlav:probSQtoQQ",       0.75  );
    parm("StringFlav:probQQ1toQQ0",     0.025 );
    parm("StringFlav:mesonUDvector",    0.5   );
    parm("StringFlav:mesonSvector",     0.6   );
    parm("StringFlav:mesonCvector",     1.5   );
    parm("StringFlav:mesonBvector",     2.5   );
    parm("StringFlav:etaSup",           0.60  );
    parm("StringFlav:etaPrimeSup",      0.15  );
    parm("StringFlav:popcornSpair",     1.0   );
    parm("StringFlav:popcornSmeson",    1.0   );
    flag("StringFlav:suppressLeadingB", false );   // kept fixed
    parm("StringZ:aLund",               0.76  );
    parm("StringZ:bLund",               0.58  );   // kept fixed
    parm("StringZ:aExtraSquark",        0.00  );   // kept fixed
    parm("StringZ:aExtraDiquark",       0.50  );   // kept fixed
    parm("StringZ:rFactC",              1.00  );   // kept fixed
    parm("StringZ:rFactB",              1.00  );   // kept fixed
    parm("StringPT:sigma",              0.36  );   // kept fixed
    parm("StringPT:enhancedFraction",   0.01  );   // kept fixed
    parm("StringPT:enhancedWidth",      2.0   );   // kept fixed
    parm("TimeShower:alphaSvalue",      0.137 );   // kept fixed
    mode("TimeShower:alphaSorder",      1     );   // kept fixed
    flag("TimeShower:alphaSuseCMW",     false );   // kept fixed
    parm("TimeShower:pTmin",            0.5   );   // kept fixed
    parm("TimeShower:pTminChgQ",        0.5   );   // kept fixed
  }

  // Full e+e- tune of flavours and FSR to LEP1 data within the
  // Rivet + Professor framework, by Hendrik Hoeth (June 2009).
  else if (eeTune == 3) {
    parm("StringFlav:probStoUD",        0.19  );
    parm("StringFlav:probQQtoQ",        0.09  );
    parm("StringFlav:probSQtoQQ",       1.00  );
    parm("StringFlav:probQQ1toQQ0",     0.027 );
    parm("StringFlav:mesonUDvector",    0.62  );
    parm("StringFlav:mesonSvector",     0.725 );
    parm("StringFlav:mesonCvector",     1.06  );
    parm("StringFlav:mesonBvector",     3.0   );
    parm("StringFlav:etaSup",           0.63  );
    parm("StringFlav:etaPrimeSup",      0.12  );
    parm("StringFlav:popcornSpair",     0.5   );   // kept fixed
    parm("StringFlav:popcornSmeson",    0.5   );   // kept fixed
    flag("StringFlav:suppressLeadingB", false );   // kept fixed
    parm("StringZ:aLund",               0.3   );   // kept fixed
    parm("StringZ:bLund",               0.8   );
    parm("StringZ:aExtraSquark",        0.00  );   // kept fixed
    parm("StringZ:aExtraDiquark",       0.50  );   // kept fixed
    parm("StringZ:rFactC",              1.00  );   // kept fixed
    parm("StringZ:rFactB",              0.67  );
    parm("StringPT:sigma",              0.304 );
    parm("StringPT:enhancedFraction",   0.01  );   // kept fixed
    parm("StringPT:enhancedWidth",      2.0   );   // kept fixed
    parm("TimeShower:alphaSvalue",      0.1383);
    mode("TimeShower:alphaSorder",      1     );   // kept fixed
    flag("TimeShower:alphaSuseCMW",     false );   // kept fixed
    parm("TimeShower:pTmin",            0.4   );   // kept fixed (near limit)
    parm("TimeShower:pTminChgQ",        0.4   );   // kept same as pTmin
  }

  // Full e+e- tune of flavours and FSR to LEP1 data, by Peter Skands
  // (September 2013). Note use of CMW convention for shower.
  else if (eeTune == 4) {
    parm("StringFlav:probStoUD",        0.21  );
    parm("StringFlav:probQQtoQ",        0.086 );
    parm("StringFlav:probSQtoQQ",       1.00  );
    parm("StringFlav:probQQ1toQQ0",     0.031 );
    parm("StringFlav:mesonUDvector",    0.45  );
    parm("StringFlav:mesonSvector",     0.60  );
    parm("StringFlav:mesonCvector",     0.95  );
    parm("StringFlav:mesonBvector",     3.0   );   // kept fixed
    parm("StringFlav:etaSup",           0.65  );
    parm("StringFlav:etaPrimeSup",      0.08  );
    parm("StringFlav:popcornSpair",     0.5   );   // kept fixed
    parm("StringFlav:popcornSmeson",    0.5   );   // kept fixed
    flag("StringFlav:suppressLeadingB", false );   // kept fixed
    parm("StringZ:aLund",               0.55  );
    parm("StringZ:bLund",               1.08  );
    parm("StringZ:aExtraSquark",        0.00  );   // kept fixed
    parm("StringZ:aExtraDiquark",       1.00  );
    parm("StringZ:rFactC",              1.00  );   // kept fixed
    parm("StringZ:rFactB",              0.85  );
    parm("StringPT:sigma",              0.305 );
    parm("StringPT:enhancedFraction",   0.01  );   // kept fixed
    parm("StringPT:enhancedWidth",      2.0   );   // kept fixed
    parm("TimeShower:alphaSvalue",      0.127 );
    mode("TimeShower:alphaSorder",      1     );   // kept fixed
    flag("TimeShower:alphaSuseCMW",     true  );
    parm("TimeShower:pTmin",            0.4   );
    parm("TimeShower:pTminChgQ",        0.4   );   // kept same as pTmin
  }

  // First e+e- tune by Nadine Fischer, using eeTune = 3 for flavour
  // composition (September 2013).
  else if (eeTune == 5) {
    parm("StringFlav:probStoUD",        0.19  );   // kept fixed
    parm("StringFlav:probQQtoQ",        0.09  );   // kept fixed
    parm("StringFlav:probSQtoQQ",       1.00  );   // kept fixed
    parm("StringFlav:probQQ1toQQ0",     0.027 );   // kept fixed
    parm("StringFlav:mesonUDvector",    0.62  );   // kept fixed
    parm("StringFlav:mesonSvector",     0.725 );   // kept fixed
    parm("StringFlav:mesonCvector",     1.06  );   // kept fixed
    parm("StringFlav:mesonBvector",     3.0   );   // kept fixed
    parm("StringFlav:etaSup",           0.63  );   // kept fixed
    parm("StringFlav:etaPrimeSup",      0.12  );   // kept fixed
    parm("StringFlav:popcornSpair",     0.5   );   // kept fixed
    parm("StringFlav:popcornSmeson",    0.5   );   // kept fixed
    flag("StringFlav:suppressLeadingB", false );   // kept fixed
    parm("StringZ:aLund",               0.386 );
    parm("StringZ:bLund",               0.977 );
    parm("StringZ:aExtraSquark",        0.00  );   // kept fixed
    parm("StringZ:aExtraDiquark",       0.940 );
    parm("StringZ:rFactC",              1.00  );   // kept fixed
    parm("StringZ:rFactB",              0.67  );   // kept fixed
    parm("StringPT:sigma",              0.286 );
    parm("StringPT:enhancedFraction",   0.01  );   // kept fixed
    parm("StringPT:enhancedWidth",      2.0   );   // kept fixed
    parm("TimeShower:alphaSvalue",      0.139 );
    mode("TimeShower:alphaSorder",      1     );   // kept fixed
    flag("TimeShower:alphaSuseCMW",     false );   // kept fixed
    parm("TimeShower:pTmin",            0.409 );
    parm("TimeShower:pTminChgQ",        0.409 );   // kept same as pTmin
  }

  // Second e+e- tune by Nadine Fischer, using eeTune = 3 for flavour
  // composition (September 2013).
  else if (eeTune == 6) {
    parm("StringFlav:probStoUD",        0.19  );   // kept fixed
    parm("StringFlav:probQQtoQ",        0.09  );   // kept fixed
    parm("StringFlav:probSQtoQQ",       1.00  );   // kept fixed
    parm("StringFlav:probQQ1toQQ0",     0.027 );   // kept fixed
    parm("StringFlav:mesonUDvector",    0.62  );   // kept fixed
    parm("StringFlav:mesonSvector",     0.725 );   // kept fixed
    parm("StringFlav:mesonCvector",     1.06  );   // kept fixed
    parm("StringFlav:mesonBvector",     3.0   );   // kept fixed
    parm("StringFlav:etaSup",           0.63  );   // kept fixed
    parm("StringFlav:etaPrimeSup",      0.12  );   // kept fixed
    parm("StringFlav:popcornSpair",     0.5   );   // kept fixed
    parm("StringFlav:popcornSmeson",    0.5   );   // kept fixed
    flag("StringFlav:suppressLeadingB", false );   // kept fixed
    parm("StringZ:aLund",               0.351 );
    parm("StringZ:bLund",               0.942 );
    parm("StringZ:aExtraSquark",        0.00  );   // kept fixed
    parm("StringZ:aExtraDiquark",       0.547 );
    parm("StringZ:rFactC",              1.00  );   // kept fixed
    parm("StringZ:rFactB",              0.67  );   // kept fixed
    parm("StringPT:sigma",              0.283 );
    parm("StringPT:enhancedFraction",   0.01  );   // kept fixed
    parm("StringPT:enhancedWidth",      2.0   );   // kept fixed
    parm("TimeShower:alphaSvalue",      0.139);
    mode("TimeShower:alphaSorder",      1     );   // kept fixed
    flag("TimeShower:alphaSuseCMW",     false );   // kept fixed
    parm("TimeShower:pTmin",            0.406 );
    parm("TimeShower:pTminChgQ",        0.406 );   // kept same as pTmin
  }

  // The Monash 2013 tune by Peter Skands, the e+e- part (January 2014).
  else if (eeTune == 7) {
    parm("StringFlav:probStoUD",        0.217 );
    parm("StringFlav:probQQtoQ",        0.081 );
    parm("StringFlav:probSQtoQQ",       0.915 );
    parm("StringFlav:probQQ1toQQ0",     0.0275);
    parm("StringFlav:mesonUDvector",    0.50  );
    parm("StringFlav:mesonSvector",     0.55  );
    parm("StringFlav:mesonCvector",     0.88  );
    parm("StringFlav:mesonBvector",     2.20  );
    parm("StringFlav:etaSup",           0.60  );
    parm("StringFlav:etaPrimeSup",      0.12  );
    parm("StringFlav:popcornSpair",     0.90  );
    parm("StringFlav:popcornSmeson",    0.50  );
    flag("StringFlav:suppressLeadingB", false );   // kept fixed
    parm("StringZ:aLund",               0.68  );
    parm("StringZ:bLund",               0.98  );
    parm("StringZ:aExtraSquark",        0.00  );   // kept fixed
    parm("StringZ:aExtraDiquark",       0.97  );
    parm("StringZ:rFactC",              1.32  );
    parm("StringZ:rFactB",              0.855 );
    parm("StringPT:sigma",              0.335 );
    parm("StringPT:enhancedFraction",   0.01  );   // kept fixed
    parm("StringPT:enhancedWidth",      2.0   );   // kept fixed
    parm("TimeShower:alphaSvalue",      0.1365);
    mode("TimeShower:alphaSorder",      1     );   // kept fixed
    flag("TimeShower:alphaSuseCMW",     false );   // kept fixed
    parm("TimeShower:pTmin",            0.5   );   // kept fixed
    parm("TimeShower:pTminChgQ",        0.5   );   // kept fixed
  }

}

//--------------------------------------------------------------------------

// Set the values related to a tune of pp/ppbar data,
// i.e. mainly for initial-state radiation and multiparton interactions.

void Settings::initTunePP( int ppTune) {

  // Restore all pp/ppbar settings to their original values.
  // Is first step for setting up a specific tune.
  if (ppTune != 0) resetTunePP();

  // Set up e+e- tune that goes with the corresponding pp tune.
  if (ppTune > 0) {
    int eeTune = 3;
    if (ppTune == 14 || ppTune >= 18) eeTune = 7;
    // The mode setting is for documentation, the real action is by initTuneEE.
    mode("Tune:ee",                            eeTune );
    initTuneEE( eeTune);
  }

  // Decide whether to use LHAPFD where possible.
  int preferLHAPDF = mode("Tune:preferLHAPDF");

  // Old ISR and MPI defaults from early and primitive comparisons with data.
  if (ppTune == 1) {
    word("PDF:pSet",                            "2"   );
    parm("SigmaProcess:alphaSvalue",            0.1265);
    flag("SigmaTotal:zeroAXB",                  true  );
    flag("SigmaDiffractive:dampen",             false );
    parm("Diffraction:largeMassSuppress",       2.0   );
    flag("TimeShower:dampenBeamRecoil",         false );
    flag("TimeShower:phiPolAsym",               false );
    parm("SpaceShower:alphaSvalue",             0.127 );
    mode("SpaceShower:alphaSorder",             1     );
    flag("SpaceShower:alphaSuseCMW",            false );
    flag("SpaceShower:samePTasMPI",             true  );
    parm("SpaceShower:pT0Ref",                  2.2   );
    parm("SpaceShower:ecmRef",                  1800.0);
    parm("SpaceShower:ecmPow",                  0.16  );
    parm("SpaceShower:pTmaxFudge",              1.0   );
    parm("SpaceShower:pTdampFudge",             1.0   );
    flag("SpaceShower:rapidityOrder",           false );
    flag("SpaceShower:phiPolAsym",              false );
    flag("SpaceShower:phiIntAsym",              false );
    parm("MultipartonInteractions:alphaSvalue", 0.127 );
    parm("MultipartonInteractions:pT0Ref",      2.15  );
    parm("MultipartonInteractions:ecmRef",      1800. );
    parm("MultipartonInteractions:ecmPow",      0.16  );
    mode("MultipartonInteractions:bProfile",    2     );
    parm("MultipartonInteractions:expPow",      1.0  );
    parm("MultipartonInteractions:a1",          0.15  );
    parm("BeamRemnants:primordialKTsoft",       0.4   );
    parm("BeamRemnants:primordialKThard",       2.1   );
    parm("BeamRemnants:halfScaleForKT",         7.0   );
    parm("BeamRemnants:halfMassForKT",          2.0   );
    mode("ColourReconnection:mode",             0     );
    parm("ColourReconnection:range",            2.5   );
  }

  // "Tune 1" simple first tune by Peter Skands to ISR and MPI, July 2009.
  else if (ppTune == 2) {
    word("PDF:pSet",                            "2"   );
    parm("SigmaProcess:alphaSvalue",            0.1265);
    flag("SigmaTotal:zeroAXB",                  true  );
    flag("SigmaDiffractive:dampen",             false );
    parm("Diffraction:largeMassSuppress",       2.0   );
    flag("TimeShower:dampenBeamRecoil",         false );
    flag("TimeShower:phiPolAsym",               false );
    parm("SpaceShower:alphaSvalue",             0.137 );
    mode("SpaceShower:alphaSorder",             1     );
    flag("SpaceShower:alphaSuseCMW",            false );
    flag("SpaceShower:samePTasMPI",             false );
    parm("SpaceShower:pT0Ref",                  2.0   );
    parm("SpaceShower:ecmRef",                  1800.0);
    parm("SpaceShower:ecmPow",                  0.0   );
    parm("SpaceShower:pTmaxFudge",              1.0   );
    parm("SpaceShower:pTdampFudge",             1.0   );
    flag("SpaceShower:rapidityOrder",           false );
    flag("SpaceShower:phiPolAsym",              false );
    flag("SpaceShower:phiIntAsym",              false );
    parm("MultipartonInteractions:alphaSvalue", 0.127 );
    parm("MultipartonInteractions:pT0Ref",      2.25  );
    parm("MultipartonInteractions:ecmRef",      1800. );
    parm("MultipartonInteractions:ecmPow",      0.24  );
    mode("MultipartonInteractions:bProfile",    1     );
    parm("MultipartonInteractions:expPow",      1.0  );
    parm("MultipartonInteractions:a1",          0.15  );
    parm("BeamRemnants:primordialKTsoft",       0.5   );
    parm("BeamRemnants:primordialKThard",       2.0   );
    parm("BeamRemnants:halfScaleForKT",         1.0   );
    parm("BeamRemnants:halfMassForKT",          1.0   );
    mode("ColourReconnection:mode",             0     );
    parm("ColourReconnection:range",            10.0  );
  }

  // Tune 2C, July 2010.
  else if (ppTune == 3) {
    word("PDF:pSet",                            "8"   );
    parm("SigmaProcess:alphaSvalue",            0.135 );
    flag("SigmaTotal:zeroAXB",                  true  );
    flag("SigmaDiffractive:dampen",             false );
    parm("Diffraction:largeMassSuppress",       2.0   );
    flag("TimeShower:dampenBeamRecoil",         true  );
    flag("TimeShower:phiPolAsym",               true  );
    parm("SpaceShower:alphaSvalue",             0.137 );
    mode("SpaceShower:alphaSorder",             1     );
    flag("SpaceShower:alphaSuseCMW",            false );
    flag("SpaceShower:samePTasMPI",             false );
    parm("SpaceShower:pT0Ref",                  2.0   );
    parm("SpaceShower:ecmRef",                  1800.0);
    parm("SpaceShower:ecmPow",                  0.0   );
    parm("SpaceShower:pTmaxFudge",              1.0   );
    parm("SpaceShower:pTdampFudge",             1.0   );
    flag("SpaceShower:rapidityOrder",           true  );
    flag("SpaceShower:phiPolAsym",              true  );
    flag("SpaceShower:phiIntAsym",              true  );
    parm("MultipartonInteractions:alphaSvalue", 0.135 );
    parm("MultipartonInteractions:pT0Ref",      2.32  );
    parm("MultipartonInteractions:ecmRef",      1800. );
    parm("MultipartonInteractions:ecmPow",      0.21  );
    mode("MultipartonInteractions:bProfile",    3     );
    parm("MultipartonInteractions:expPow",      1.6   );
    parm("MultipartonInteractions:a1",          0.15  );
    parm("BeamRemnants:primordialKTsoft",       0.5   );
    parm("BeamRemnants:primordialKThard",       2.0   );
    parm("BeamRemnants:halfScaleForKT",         1.0   );
    parm("BeamRemnants:halfMassForKT",          1.0   );
    mode("ColourReconnection:mode",             0     );
    parm("ColourReconnection:range",            3.0   );
  }

  // Tune 2M, July 2010.
  else if (ppTune == 4) {
    word("PDF:pSet",                            "4"   );
    parm("SigmaProcess:alphaSvalue",            0.1265);
    flag("SigmaTotal:zeroAXB",                  true  );
    flag("SigmaDiffractive:dampen",             false );
    parm("Diffraction:largeMassSuppress",       2.0   );
    flag("TimeShower:dampenBeamRecoil",         true  );
    flag("TimeShower:phiPolAsym",               true  );
    parm("SpaceShower:alphaSvalue",             0.130 );
    mode("SpaceShower:alphaSorder",             1     );
    flag("SpaceShower:alphaSuseCMW",            false );
    flag("SpaceShower:samePTasMPI",             false );
    parm("SpaceShower:pT0Ref",                  2.0   );
    parm("SpaceShower:ecmRef",                  1800.0);
    parm("SpaceShower:ecmPow",                  0.0   );
    parm("SpaceShower:pTmaxFudge",              1.0   );
    parm("SpaceShower:pTdampFudge",             1.0   );
    flag("SpaceShower:rapidityOrder",           true  );
    flag("SpaceShower:phiPolAsym",              true  );
    flag("SpaceShower:phiIntAsym",              true  );
    parm("MultipartonInteractions:alphaSvalue", 0.127 );
    parm("MultipartonInteractions:pT0Ref",      2.455 );
    parm("MultipartonInteractions:ecmRef",      1800. );
    parm("MultipartonInteractions:ecmPow",      0.26  );
    mode("MultipartonInteractions:bProfile",    3     );
    parm("MultipartonInteractions:expPow",      1.15  );
    parm("MultipartonInteractions:a1",          0.15  );
    parm("BeamRemnants:primordialKTsoft",       0.5   );
    parm("BeamRemnants:primordialKThard",       2.0   );
    parm("BeamRemnants:halfScaleForKT",         1.0   );
    parm("BeamRemnants:halfMassForKT",          1.0   );
    mode("ColourReconnection:mode",             0     );
    parm("ColourReconnection:range",            3.0   );
  }

  // Tune 4C, October 2010.
  else if (ppTune == 5) {
    word("PDF:pSet",                            "8"   );
    parm("SigmaProcess:alphaSvalue",            0.135 );
    flag("SigmaTotal:zeroAXB",                  true  );
    flag("SigmaDiffractive:dampen",             true  );
    parm("SigmaDiffractive:maxXB",              65.0  );
    parm("SigmaDiffractive:maxAX",              65.0  );
    parm("SigmaDiffractive:maxXX",              65.0  );
    parm("Diffraction:largeMassSuppress",       2.0   );
    flag("TimeShower:dampenBeamRecoil",         true  );
    flag("TimeShower:phiPolAsym",               true  );
    parm("SpaceShower:alphaSvalue",             0.137 );
    mode("SpaceShower:alphaSorder",             1     );
    flag("SpaceShower:alphaSuseCMW",            false );
    flag("SpaceShower:samePTasMPI",             false );
    parm("SpaceShower:pT0Ref",                  2.0   );
    parm("SpaceShower:ecmRef",                  1800.0);
    parm("SpaceShower:ecmPow",                  0.0   );
    parm("SpaceShower:pTmaxFudge",              1.0   );
    parm("SpaceShower:pTdampFudge",             1.0   );
    flag("SpaceShower:rapidityOrder",           true  );
    flag("SpaceShower:phiPolAsym",              true  );
    flag("SpaceShower:phiIntAsym",              true  );
    parm("MultipartonInteractions:alphaSvalue", 0.135 );
    parm("MultipartonInteractions:pT0Ref",      2.085 );
    parm("MultipartonInteractions:ecmRef",      1800. );
    parm("MultipartonInteractions:ecmPow",      0.19  );
    mode("MultipartonInteractions:bProfile",    3     );
    parm("MultipartonInteractions:expPow",      2.0   );
    parm("MultipartonInteractions:a1",          0.15  );
    parm("BeamRemnants:primordialKTsoft",       0.5   );
    parm("BeamRemnants:primordialKThard",       2.0   );
    parm("BeamRemnants:halfScaleForKT",         1.0   );
    parm("BeamRemnants:halfMassForKT",          1.0   );
    mode("ColourReconnection:mode",             0     );
    parm("ColourReconnection:range",            1.5   );
  }

  // Tune 4Cx, January 2011.
  else if (ppTune == 6) {
    word("PDF:pSet",                            "8"   );
    parm("SigmaProcess:alphaSvalue",            0.135 );
    flag("SigmaTotal:zeroAXB",                  true  );
    flag("SigmaDiffractive:dampen",             true  );
    parm("SigmaDiffractive:maxXB",              65.0  );
    parm("SigmaDiffractive:maxAX",              65.0  );
    parm("SigmaDiffractive:maxXX",              65.0  );
    parm("Diffraction:largeMassSuppress",       2.0   );
    flag("TimeShower:dampenBeamRecoil",         true  );
    flag("TimeShower:phiPolAsym",               true  );
    parm("SpaceShower:alphaSvalue",             0.137 );
    mode("SpaceShower:alphaSorder",             1     );
    flag("SpaceShower:alphaSuseCMW",            false );
    flag("SpaceShower:samePTasMPI",             false );
    parm("SpaceShower:pT0Ref",                  2.0   );
    parm("SpaceShower:ecmRef",                  1800.0);
    parm("SpaceShower:ecmPow",                  0.0   );
    parm("SpaceShower:pTmaxFudge",              1.0   );
    parm("SpaceShower:pTdampFudge",             1.0   );
    flag("SpaceShower:rapidityOrder",           true  );
    flag("SpaceShower:phiPolAsym",              true  );
    flag("SpaceShower:phiIntAsym",              true  );
    parm("MultipartonInteractions:alphaSvalue", 0.135 );
    parm("MultipartonInteractions:pT0Ref",      2.15  );
    parm("MultipartonInteractions:ecmRef",      1800. );
    parm("MultipartonInteractions:ecmPow",      0.19  );
    mode("MultipartonInteractions:bProfile",    4     );
    parm("MultipartonInteractions:expPow",      1.0   );
    parm("MultipartonInteractions:a1",          0.15  );
    parm("BeamRemnants:primordialKTsoft",       0.5   );
    parm("BeamRemnants:primordialKThard",       2.0   );
    parm("BeamRemnants:halfScaleForKT",         1.0   );
    parm("BeamRemnants:halfMassForKT",          1.0   );
    mode("ColourReconnection:mode",             0     );
    parm("ColourReconnection:range",            1.5   );
  }

  // The Monash 2013 tune by Peter Skands, the pp part (January 2014).
  else if (ppTune == 14) {
    word("PDF:pSet",                            "13"  );   // NNPDF
    parm("SigmaProcess:alphaSvalue",            0.130 );   // same as PDF
    flag("SigmaTotal:zeroAXB",                  true  );
    flag("SigmaDiffractive:dampen",             true  );
    parm("SigmaDiffractive:maxXB",              65.0  );
    parm("SigmaDiffractive:maxAX",              65.0  );
    parm("SigmaDiffractive:maxXX",              65.0  );
    parm("Diffraction:largeMassSuppress",       4.0   );
    flag("TimeShower:dampenBeamRecoil",         true  );
    flag("TimeShower:phiPolAsym",               true  );
    parm("SpaceShower:alphaSvalue",             0.1365);   // same as FSR
    mode("SpaceShower:alphaSorder",             1     );
    flag("SpaceShower:alphaSuseCMW",            false );
    flag("SpaceShower:samePTasMPI",             false );
    parm("SpaceShower:pT0Ref",                  2.0   );
    parm("SpaceShower:ecmRef",                  7000.0);
    parm("SpaceShower:ecmPow",                  0.0   );
    parm("SpaceShower:pTmaxFudge",              1.0   );
    parm("SpaceShower:pTdampFudge",             1.0   );
    flag("SpaceShower:rapidityOrder",           true  );
    flag("SpaceShower:phiPolAsym",              true  );
    flag("SpaceShower:phiIntAsym",              true  );
    parm("MultipartonInteractions:alphaSvalue", 0.130 );   // same as PDF
    parm("MultipartonInteractions:pT0Ref",      2.28  );
    parm("MultipartonInteractions:ecmRef",      7000. );
    parm("MultipartonInteractions:ecmPow",      0.215 );
    mode("MultipartonInteractions:bProfile",    3     );
    parm("MultipartonInteractions:expPow",      1.85  );
    parm("MultipartonInteractions:a1",          0.15  );
    parm("BeamRemnants:primordialKTsoft",       0.9   );
    parm("BeamRemnants:primordialKThard",       1.8   );
    parm("BeamRemnants:halfScaleForKT",         1.5   );
    parm("BeamRemnants:halfMassForKT",          1.0   );
    mode("ColourReconnection:mode",             0     );
    parm("ColourReconnection:range",            1.80  );
  }

  // Several ATLAS and CMS tunes start out from Tune 4C.
  else if (ppTune < 18) {
    parm("SigmaProcess:alphaSvalue",            0.135 );
    flag("SigmaTotal:zeroAXB",                  true  );
    flag("SigmaDiffractive:dampen",             true  );
    parm("SigmaDiffractive:maxXB",              65.0  );
    parm("SigmaDiffractive:maxAX",              65.0  );
    parm("SigmaDiffractive:maxXX",              65.0  );
    parm("Diffraction:largeMassSuppress",       2.0   );
    flag("TimeShower:dampenBeamRecoil",         true  );
    flag("TimeShower:phiPolAsym",               true  );
    parm("SpaceShower:alphaSvalue",             0.137 );
    mode("SpaceShower:alphaSorder",             1     );
    flag("SpaceShower:alphaSuseCMW",            false );
    flag("SpaceShower:samePTasMPI",             false );
    parm("SpaceShower:pT0Ref",                  2.0   );
    parm("SpaceShower:ecmRef",                  1800.0);
    parm("SpaceShower:ecmPow",                  0.0   );
    parm("SpaceShower:pTmaxFudge",              1.0   );
    parm("SpaceShower:pTdampFudge",             1.0   );
    flag("SpaceShower:rapidityOrder",           true );
    flag("SpaceShower:phiPolAsym",              true  );
    flag("SpaceShower:phiIntAsym",              true  );
    parm("MultipartonInteractions:alphaSvalue", 0.135 );
    parm("MultipartonInteractions:pT0Ref",      2.085 );
    parm("MultipartonInteractions:ecmRef",      1800. );
    parm("MultipartonInteractions:ecmPow",      0.19  );
    mode("MultipartonInteractions:bProfile",    3     );
    parm("MultipartonInteractions:expPow",      2.0   );
    parm("MultipartonInteractions:a1",          0.15  );
    parm("BeamRemnants:primordialKTsoft",       0.5   );
    parm("BeamRemnants:primordialKThard",       2.0   );
    parm("BeamRemnants:halfScaleForKT",         1.0   );
    parm("BeamRemnants:halfMassForKT",          1.0   );
    mode("ColourReconnection:mode",             0     );
    parm("ColourReconnection:range",            1.5   );

    // Several ATLAS tunes in the A2 and AU2 series, see
    // ATLAS note ATL-PHYS-PUB-2012-003 (August 2012).
    // ATLAS MB tune A2-CTEQ6L1.
    if (ppTune == 7) {
      if (preferLHAPDF == 1)
        word("PDF:pSet",       "LHAPDF5:cteq6ll.LHpdf");
      else if (preferLHAPDF == 2)
        word("PDF:pSet",             "LHAPDF6:cteq6l1");
      else word("PDF:pSet",                     "8"   );
      flag("SpaceShower:rapidityOrder",         false );
      parm("MultipartonInteractions:pT0Ref",    2.18  );
      parm("MultipartonInteractions:ecmPow",    0.22  );
      mode("MultipartonInteractions:bProfile",  4     );
      parm("MultipartonInteractions:expPow",    1.0   );
      parm("MultipartonInteractions:a1",        0.06  );
      parm("ColourReconnection:range",          1.55  );
    }

    // ATLAS MB tune A2-MSTW2008LO.
    else if (ppTune == 8) {
      if (preferLHAPDF == 1)
        word("PDF:pSet", "LHAPDF5:MSTW2008lo68cl.LHgrid");
      else if (preferLHAPDF == 2)
        word("PDF:pSet",      "LHAPDF6:MSTW2008lo68cl");
      else word("PDF:pSet",                     "5"   );
      flag("SpaceShower:rapidityOrder",         false );
      parm("MultipartonInteractions:pT0Ref",    1.90  );
      parm("MultipartonInteractions:ecmPow",    0.30  );
      mode("MultipartonInteractions:bProfile",  4     );
      parm("MultipartonInteractions:expPow",    1.0   );
      parm("MultipartonInteractions:a1",        0.03  );
      parm("ColourReconnection:range",          2.28  );
    }

    // ATLAS UE tune AU2-CTEQ6L1.
    if (ppTune == 9) {
      if (preferLHAPDF == 1)
        word("PDF:pSet",       "LHAPDF5:cteq6ll.LHpdf");
      else if (preferLHAPDF == 2)
        word("PDF:pSet",             "LHAPDF6:cteq6l1");
      else word("PDF:pSet",                     "8"   );
      flag("SpaceShower:rapidityOrder",         false );
      parm("MultipartonInteractions:pT0Ref",    2.13  );
      parm("MultipartonInteractions:ecmPow",    0.21  );
      mode("MultipartonInteractions:bProfile",  4     );
      parm("MultipartonInteractions:expPow",    1.0   );
      parm("MultipartonInteractions:a1",        0.00  );
      parm("ColourReconnection:range",          2.21  );
    }

    // ATLAS UE tune AU2-MSTW2008LO.
    else if (ppTune == 10) {
      if (preferLHAPDF == 1)
        word("PDF:pSet", "LHAPDF5:MSTW2008lo68cl.LHgrid");
      else if (preferLHAPDF == 2)
        word("PDF:pSet",      "LHAPDF6:MSTW2008lo68cl");
      else word("PDF:pSet",                     "5"   );
      flag("SpaceShower:rapidityOrder",         false );
      parm("MultipartonInteractions:pT0Ref",    1.87  );
      parm("MultipartonInteractions:ecmPow",    0.28  );
      mode("MultipartonInteractions:bProfile",  4     );
      parm("MultipartonInteractions:expPow",    1.0   );
      parm("MultipartonInteractions:a1",        0.01  );
      parm("ColourReconnection:range",          5.32  );
    }

    // ATLAS UE tune AU2-CT10.
    else if (ppTune == 11) {
      if (preferLHAPDF == 2)
        word("PDF:pSet",                "LHAPDF6:CT10");
      else
        word("PDF:pSet",         "LHAPDF5:CT10.LHgrid");
      flag("SpaceShower:rapidityOrder",         false );
      parm("MultipartonInteractions:pT0Ref",    1.70  );
      parm("MultipartonInteractions:ecmPow",    0.16  );
      mode("MultipartonInteractions:bProfile",  4     );
      parm("MultipartonInteractions:expPow",    1.0   );
      parm("MultipartonInteractions:a1",        0.10  );
      parm("ColourReconnection:range",          4.67  );
    }

    // ATLAS UE tune AU2-MRST2007LO*.
    else if (ppTune == 12) {
      if (preferLHAPDF == 1)
        word("PDF:pSet", "LHAPDF5:MRST2007lomod.LHgrid");
      else if (preferLHAPDF == 2)
        word("PDF:pSet",       "LHAPDF6:MRST2007lomod");
      else word("PDF:pSet",                     "3"   );
      flag("SpaceShower:rapidityOrder",         false );
      parm("MultipartonInteractions:pT0Ref",    2.39  );
      parm("MultipartonInteractions:ecmPow",    0.24  );
      mode("MultipartonInteractions:bProfile",  4     );
      parm("MultipartonInteractions:expPow",    1.0   );
      parm("MultipartonInteractions:a1",        0.01  );
      parm("ColourReconnection:range",          1.76  );
    }

    // ATLAS UE tune AU2-MRST2007LO**.
    else if (ppTune == 13) {
      if (preferLHAPDF == 1)
        word("PDF:pSet",     "LHAPDF5:MRSTMCal.LHgrid");
      else if (preferLHAPDF == 2)
        word("PDF:pSet",            "LHAPDF6:MRSTMCal");
      else word("PDF:pSet",                     "4"   );
      flag("SpaceShower:rapidityOrder",         false );
      parm("MultipartonInteractions:pT0Ref",    2.57  );
      parm("MultipartonInteractions:ecmPow",    0.23  );
      mode("MultipartonInteractions:bProfile",  4     );
      parm("MultipartonInteractions:expPow",    1.0   );
      parm("MultipartonInteractions:a1",        0.01  );
      parm("ColourReconnection:range",          1.47  );
    }

    // The CMS UE tunes CUETP8S1-CTEQ6L1 and CUETP8S1-HERAPDF1.5LO,
    // see the note CMS PAS GEN-14-001 (April 2014).
    // CMS UE tune CUETP8S1-CTEQ6L1.
    else if (ppTune == 15) {
      if (preferLHAPDF == 1)
        word("PDF:pSet",       "LHAPDF5:cteq6ll.LHpdf");
      else if (preferLHAPDF == 2)
        word("PDF:pSet",             "LHAPDF6:cteq6l1");
      else word("PDF:pSet",                     "8"   );
      parm("MultipartonInteractions:pT0Ref",    2.1006);
      parm("MultipartonInteractions:ecmPow",    0.2106);
      parm("MultipartonInteractions:expPow",    1.6089);
      parm("MultipartonInteractions:a1",        0.00  );
      parm("ColourReconnection:range",          3.3126);
    }

    // CMS UE tune CUETP8S1-HERAPDF1.5LO.
    else if (ppTune == 16) {
      if (preferLHAPDF == 2)
        word("PDF:pSet",     "LHAPDF6:HERAPDF15LO_EIG");
      else
        word("PDF:pSet", "LHAPDF5:HERAPDF1.5LO_EIG.LHgrid");
      parm("MultipartonInteractions:pT0Ref",    2.0001);
      parm("MultipartonInteractions:ecmPow",    0.2499);
      parm("MultipartonInteractions:expPow",    1.6905);
      parm("MultipartonInteractions:a1",        0.00  );
      parm("ColourReconnection:range",          6.0964);
    }

    // ATLAS tune AZ to the Z0/gamma* pTspectrum, see the note
    // CERN-PH-EP-2014-075 [arXiv:1406.3660 [hep-ex]] (June 2014).
    else if (ppTune == 17) {
      parm("SpaceShower:alphaSvalue",           0.1237);
      parm("SpaceShower:pT0Ref",                0.59  );
      parm("MultipartonInteractions:pT0Ref",    2.18  );
      parm("BeamRemnants:primordialKThard",     2.0   );
    }
  }

  // Several ATLAS and CMS tunes start out from Monash 2013 tune.
  else if (ppTune >= 18) {
    word("PDF:pSet",                            "13"  );   // NNPDF
    parm("SigmaProcess:alphaSvalue",            0.130 );   // same as PDF
    flag("SigmaTotal:zeroAXB",                  true  );
    flag("SigmaDiffractive:dampen",             true  );
    parm("SigmaDiffractive:maxXB",              65.0  );
    parm("SigmaDiffractive:maxAX",              65.0  );
    parm("SigmaDiffractive:maxXX",              65.0  );
    parm("Diffraction:largeMassSuppress",       4.0   );
    flag("TimeShower:dampenBeamRecoil",         true  );
    flag("TimeShower:phiPolAsym",               true  );
    parm("SpaceShower:alphaSvalue",             0.1365);   // same as FSR
    mode("SpaceShower:alphaSorder",             1     );
    flag("SpaceShower:alphaSuseCMW",            false );
    flag("SpaceShower:samePTasMPI",             false );
    parm("SpaceShower:pT0Ref",                  2.0   );
    parm("SpaceShower:ecmRef",                  7000.0);
    parm("SpaceShower:ecmPow",                  0.0   );
    parm("SpaceShower:pTmaxFudge",              1.0   );
    parm("SpaceShower:pTdampFudge",             1.0   );
    flag("SpaceShower:rapidityOrder",           true  );
    flag("SpaceShower:phiPolAsym",              true  );
    flag("SpaceShower:phiIntAsym",              true  );
    parm("MultipartonInteractions:alphaSvalue", 0.130 );   // same as PDF
    parm("MultipartonInteractions:pT0Ref",      2.28  );
    parm("MultipartonInteractions:ecmRef",      7000. );
    parm("MultipartonInteractions:ecmPow",      0.215 );
    mode("MultipartonInteractions:bProfile",    3     );
    parm("MultipartonInteractions:expPow",      1.85  );
    parm("MultipartonInteractions:a1",          0.15  );
    parm("BeamRemnants:primordialKTsoft",       0.9   );
    parm("BeamRemnants:primordialKThard",       1.8   );
    parm("BeamRemnants:halfScaleForKT",         1.5   );
    parm("BeamRemnants:halfMassForKT",          1.0   );
    mode("ColourReconnection:mode",             0     );
    parm("ColourReconnection:range",            1.80  );

    // CMS tune MonashStar = CUETP8M1-NNPDF2.3LO.
    // See R.D. Field, presentation at MPI@LHC 2014, Krakow, Poland.
    if (ppTune == 18) {
      parm("MultipartonInteractions:pT0Ref",    2.4024);
      parm("MultipartonInteractions:ecmPow",    0.2521);
      parm("MultipartonInteractions:expPow",    1.60  );
    }

    // The ATLAS A14 tunes, central tune with CTEQL1.
    // See ATL-PHYS-PUB-2014-021 (November 2014).
    // Warning: note that TimeShower:alphaSvalue is set here, although
    // normally it would be in the domain of ee tunes. This makes the
    // order of Tune:ee and Tune:pp commands relevant.
    else if (ppTune == 19) {
      if (preferLHAPDF == 1)
        word("PDF:pSet",       "LHAPDF5:cteq6ll.LHpdf");
      else if (preferLHAPDF == 2)
        word("PDF:pSet",             "LHAPDF6:cteq6l1");
      else word("PDF:pSet",                     "8"   );
      parm("SigmaProcess:alphaSvalue",          0.144 );
      parm("TimeShower:alphaSvalue",            0.126 );
      parm("SpaceShower:alphaSvalue",           0.125 );
      parm("SpaceShower:pT0Ref",                1.3   );
      parm("SpaceShower:pTmaxFudge",            0.95   );
      parm("SpaceShower:pTdampFudge",           1.21  );
      parm("MultipartonInteractions:alphaSvalue",0.118);
      parm("MultipartonInteractions:pT0Ref",    1.98  );
      parm("BeamRemnants:primordialKThard",     1.72  );
      parm("ColourReconnection:range",          2.08  );
    }

    // The ATLAS A14 tunes, central tune with MSTW2008LO.
    else if (ppTune == 20) {
      if (preferLHAPDF == 1)
        word("PDF:pSet", "LHAPDF5:MSTW2008lo68cl.LHgrid");
      else if (preferLHAPDF == 2)
        word("PDF:pSet",      "LHAPDF6:MSTW2008lo68cl");
      else word("PDF:pSet",                     "5"   );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.129 );
      parm("SpaceShower:alphaSvalue",           0.129 );
      parm("SpaceShower:pT0Ref",                1.62  );
      parm("SpaceShower:pTmaxFudge",            0.92  );
      parm("SpaceShower:pTdampFudge",           1.14  );
      parm("MultipartonInteractions:alphaSvalue",0.130);
      parm("MultipartonInteractions:pT0Ref",    2.28  );
      parm("BeamRemnants:primordialKThard",     1.82  );
      parm("ColourReconnection:range",          1.87  );
    }

    // The ATLAS A14 tunes, central tune with NNPDF2.3LO.
    else if (ppTune == 21) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.127 );
      parm("SpaceShower:alphaSvalue",           0.127 );
      parm("SpaceShower:pT0Ref",                1.56  );
      parm("SpaceShower:pTmaxFudge",            0.91  );
      parm("SpaceShower:pTdampFudge",           1.05  );
      parm("MultipartonInteractions:alphaSvalue",0.126);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.71  );
    }

    // The ATLAS A14 tunes, central tune with HERAPDF1.5LO.
    else if (ppTune == 22) {
      if (preferLHAPDF == 2)
        word("PDF:pSet",     "LHAPDF6:HERAPDF15LO_EIG");
      else
        word("PDF:pSet", "LHAPDF5:HERAPDF1.5LO_EIG.LHgrid");
      parm("SigmaProcess:alphaSvalue",          0.141 );
      parm("TimeShower:alphaSvalue",            0.130 );
      parm("SpaceShower:alphaSvalue",           0.128);
      parm("SpaceShower:pT0Ref",                1.61  );
      parm("SpaceShower:pTmaxFudge",            0.95  );
      parm("SpaceShower:pTdampFudge",           1.10  );
      parm("MultipartonInteractions:alphaSvalue",0.123);
      parm("MultipartonInteractions:pT0Ref",    2.14  );
      parm("BeamRemnants:primordialKThard",     1.83  );
      parm("ColourReconnection:range",          1.78  );
    }

    // The ATLAS A14 tunes, variation 1+.
    else if (ppTune == 23) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.127 );
      parm("SpaceShower:alphaSvalue",           0.127 );
      parm("SpaceShower:pT0Ref",                1.56  );
      parm("SpaceShower:pTmaxFudge",            0.91  );
      parm("SpaceShower:pTdampFudge",           1.05  );
      parm("MultipartonInteractions:alphaSvalue",0.131);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.73  );
    }

    // The ATLAS A14 tunes, variation 1-.
    else if (ppTune == 24) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.127 );
      parm("SpaceShower:alphaSvalue",           0.127 );
      parm("SpaceShower:pT0Ref",                1.56  );
      parm("SpaceShower:pTmaxFudge",            0.91  );
      parm("SpaceShower:pTdampFudge",           1.05  );
      parm("MultipartonInteractions:alphaSvalue",0.121);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.69  );
    }

    // The ATLAS A14 tunes, variation 2+.
    else if (ppTune == 25) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.139 );
      parm("SpaceShower:alphaSvalue",           0.127 );
      parm("SpaceShower:pT0Ref",                1.60  );
      parm("SpaceShower:pTmaxFudge",            0.91  );
      parm("SpaceShower:pTdampFudge",           1.04  );
      parm("MultipartonInteractions:alphaSvalue",0.126);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.71  );
    }

    // The ATLAS A14 tunes, variation 2-.
    else if (ppTune == 26) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.111 );
      parm("SpaceShower:alphaSvalue",           0.127 );
      parm("SpaceShower:pT0Ref",                1.50  );
      parm("SpaceShower:pTmaxFudge",            0.91  );
      parm("SpaceShower:pTdampFudge",           1.08  );
      parm("MultipartonInteractions:alphaSvalue",0.126);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.71  );
    }

    // The ATLAS A14 tunes, variation 3a+.
    else if (ppTune == 27) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.136 );
      parm("SpaceShower:alphaSvalue",           0.127 );
      parm("SpaceShower:pT0Ref",                1.67  );
      parm("SpaceShower:pTmaxFudge",            0.98  );
      parm("SpaceShower:pTdampFudge",           1.36  );
      parm("MultipartonInteractions:alphaSvalue",0.125);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.71  );
    }

    // The ATLAS A14 tunes, variation 3a-.
    else if (ppTune == 28) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.124 );
      parm("SpaceShower:alphaSvalue",           0.127 );
      parm("SpaceShower:pT0Ref",                1.51  );
      parm("SpaceShower:pTmaxFudge",            0.88  );
      parm("SpaceShower:pTdampFudge",           0.93  );
      parm("MultipartonInteractions:alphaSvalue",0.127);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.71  );
    }

    // The ATLAS A14 tunes, variation 3b+.
    else if (ppTune == 29) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.114 );
      parm("SpaceShower:alphaSvalue",           0.129 );
      parm("SpaceShower:pT0Ref",                1.56  );
      parm("SpaceShower:pTmaxFudge",            1.00  );
      parm("SpaceShower:pTdampFudge",           1.04  );
      parm("MultipartonInteractions:alphaSvalue",0.126);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.71  );
    }

    // The ATLAS A14 tunes, variation 3b-.
    else if (ppTune == 30) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.138 );
      parm("SpaceShower:alphaSvalue",           0.126 );
      parm("SpaceShower:pT0Ref",                1.56  );
      parm("SpaceShower:pTmaxFudge",            0.83  );
      parm("SpaceShower:pTdampFudge",           1.07  );
      parm("MultipartonInteractions:alphaSvalue",0.126);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.71  );
    }

    // The ATLAS A14 tunes, variation 3c+.
    else if (ppTune == 31) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.127 );
      parm("SpaceShower:alphaSvalue",           0.140 );
      parm("SpaceShower:pT0Ref",                1.56  );
      parm("SpaceShower:pTmaxFudge",            0.91  );
      parm("SpaceShower:pTdampFudge",           1.05  );
      parm("MultipartonInteractions:alphaSvalue",0.126);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.71  );
    }

    // The ATLAS A14 tunes, variation 3c-.
    else if (ppTune == 32) {
      word("PDF:pSet",                          "13"  );
      parm("SigmaProcess:alphaSvalue",          0.140 );
      parm("TimeShower:alphaSvalue",            0.127 );
      parm("SpaceShower:alphaSvalue",           0.115 );
      parm("SpaceShower:pT0Ref",                1.56  );
      parm("SpaceShower:pTmaxFudge",            0.91  );
      parm("SpaceShower:pTdampFudge",           1.05  );
      parm("MultipartonInteractions:alphaSvalue",0.126);
      parm("MultipartonInteractions:pT0Ref",    2.09  );
      parm("BeamRemnants:primordialKThard",     1.88  );
      parm("ColourReconnection:range",          1.71  );
    }

  }

}

//--------------------------------------------------------------------------

// Convert string to lowercase for case-insensitive comparisons.
// Also remove initial and trailing blanks, if any.

string Settings::toLower(const string& name) {

  // Copy string without initial and trailing blanks.
  if (name.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return "";
  int firstChar = name.find_first_not_of(" \n\t\v\b\r\f\a");
  int lastChar  = name.find_last_not_of(" \n\t\v\b\r\f\a");
  string temp   = name.substr( firstChar, lastChar + 1 - firstChar);

  // Convert to lowercase letter by letter.
  for (int i = 0; i < int(temp.length()); ++i) temp[i] = tolower(temp[i]);
  return temp;

}

//--------------------------------------------------------------------------

// Allow several alternative inputs for true/false.

bool Settings::boolString(string tag) {

  string tagLow = toLower(tag);
  return ( tagLow == "true" || tagLow == "1" || tagLow == "on"
  || tagLow == "yes" || tagLow == "ok" );

}

//--------------------------------------------------------------------------

// Extract XML value string following XML attribute.

string Settings::attributeValue(string line, string attribute) {

  if (line.find(attribute) == string::npos) return "";
  int iBegAttri = line.find(attribute);
  int iBegQuote = line.find("\"", iBegAttri + 1);
  int iEndQuote = line.find("\"", iBegQuote + 1);
  return line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);

}

//--------------------------------------------------------------------------

// Extract XML bool value following XML attribute.

bool Settings::boolAttributeValue(string line, string attribute) {

  string valString = attributeValue(line, attribute);
  if (valString == "") return false;
  return boolString(valString);

}

//--------------------------------------------------------------------------

// Extract XML int value following XML attribute.

int Settings::intAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0;
  istringstream valStream(valString);
  int intVal;
  valStream >> intVal;
  return intVal;

}

//--------------------------------------------------------------------------

// Extract XML double value following XML attribute.

double Settings::doubleAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0.;
  istringstream valStream(valString);
  double doubleVal;
  valStream >> doubleVal;
  return doubleVal;

}

//--------------------------------------------------------------------------

// Extract XML bool vector value following XML attribute.

vector<bool> Settings::boolVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return vector<bool>(1, false);
  vector<bool> vectorVal;
  size_t       stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    vectorVal.push_back(boolString(valStream.str()));
  }
  return vectorVal;

}

//--------------------------------------------------------------------------

// Extract XML int vector value following XML attribute.

vector<int> Settings::intVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return vector<int>(1, 0);
  int         intVal;
  vector<int> vectorVal;
  size_t      stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    valStream >> intVal;
    vectorVal.push_back(intVal);
  }
  return vectorVal;

}

//--------------------------------------------------------------------------

// Extract XML double vector value following XML attribute.

vector<double> Settings::doubleVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return vector<double>(1, 0.);
  double         doubleVal;
  vector<double> vectorVal;
  size_t         stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    valStream >> doubleVal;
    vectorVal.push_back(doubleVal);
  }
  return vectorVal;

}

//==========================================================================

} // end namespace Pythia8
