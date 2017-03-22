// Settings.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Settings class.

#include "Settings.h"

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
         && tag != "<parm" && tag != "<parmfix" && tag != "<vect"
	 && tag != "<vectfix" && tag != "<word" 
         && tag != "<wordfix" && tag != "<aidx") continue;

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
        int value  = intAttributeValue( line, "default="); 
        int minVal = intAttributeValue( line, "min=");
        int maxVal = intAttributeValue( line, "max=");
        addMode( name, value, hasMin, hasMax, minVal, maxVal);	  
    
      // Check for occurence of a double and add to parm map.
      } else if (tag == "<parm" || tag == "<parmfix") {
        double value  = doubleAttributeValue( line, "default="); 
        double minVal = doubleAttributeValue( line, "min=");
        double maxVal = doubleAttributeValue( line, "max=");
        addParm( name, value, hasMin, hasMax, minVal, maxVal);
	
      // Check for occurence of a vector and add to vect map.
      } else if (tag == "<vect" || tag == "<vectfix") {
        vector<double> value = vectorAttributeValue( line, "default="); 
        double minVal = doubleAttributeValue( line, "min=");
        double maxVal = doubleAttributeValue( line, "max=");
        addVect( name, value, hasMin, hasMax, minVal, maxVal);

      // Check for occurence of a string and add to word map.
      } else if (tag == "<word" || tag == "<wordfix") {
        string value = attributeValue( line, "default="); 
        addWord( name, value);
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
  vects.clear();
  words.clear();

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
  else if (isVect(name)) inDataBase = 4; 
  else if (isWord(name)) inDataBase = 5; 

  // For backwards compatibility: multiple -> multiparton, MI -> MPI.
  if (inDataBase == 0) {
    bool retry = false;
    string nameLower = toLower(name);
    if (nameLower.find("multiple") != string::npos) {
      int firstMI = nameLower.find_first_of("multiple");
      name.replace(firstMI, 8, "Multiparton");  
      retry = true; 
    }
    if (!retry && nameLower.find("mi") != string::npos) {
      int firstMI = nameLower.find_first_of("mi");
      name.replace(firstMI, 2, "MPI");   
      retry = true; 
    }
    if (retry) {
      if      (isFlag(name)) inDataBase = 1;   
      else if (isMode(name)) inDataBase = 2;   
      else if (isParm(name)) inDataBase = 3; 
      else if (isVect(name)) inDataBase = 4; 
      else if (isWord(name)) inDataBase = 5; 
    }
  }

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
    mode(name, value);
        
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
        
  // Update vect map.
  } else if (inDataBase == 4) {
    istringstream vectData(valueString);
    vector<double> value(vectorAttributeValue(
      "value=\"" + valueString + "\"", "value="));
    if (!vectData) {
      if (warn) os << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;  
    }  
    vect(name, value);
        
  // Update word map.
  } else {
    word(name, valueString);
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
  map<string, Vect>::iterator vectEntry = vects.begin();
  map<string, Word>::iterator wordEntry = words.begin();

  // Loop while there is something left to do.
  while (flagEntry != flags.end() || modeEntry != modes.end() 
    || parmEntry != parms.end() || vectEntry != vects.end() 
    || wordEntry != words.end()) {

    // Check if a flag is next in lexigraphical order; if so print it.
    if ( flagEntry != flags.end() 
      && ( modeEntry == modes.end() || flagEntry->first < modeEntry->first ) 
      && ( parmEntry == parms.end() || flagEntry->first < parmEntry->first )
      && ( vectEntry == vects.end() || flagEntry->first < vectEntry->first )
      && ( wordEntry == words.end() || flagEntry->first < wordEntry->first ) 
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
      && ( vectEntry == vects.end() || modeEntry->first < vectEntry->first ) 
      && ( wordEntry == words.end() || modeEntry->first < wordEntry->first ) 
      ) {
      int valNow = modeEntry->second.valNow;
      int valDefault = modeEntry->second.valDefault;
      if ( writeAll || valNow != valDefault ) 
        os << modeEntry->second.name << " = " << valNow << "\n"; 
      ++modeEntry;
      
    // Else check if parm is next, and if so print it; 
    // fixed or scientific depending on value.
    } else if ( parmEntry != parms.end()
      && ( vectEntry == vects.end() || parmEntry->first < vectEntry->first ) 
      && ( wordEntry == words.end() || parmEntry->first < wordEntry->first ) 
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

    // Else check if vect is next, and if so print it; 
    // fixed or scientific depending on value.
    } else if ( vectEntry != vects.end()
      && ( wordEntry == words.end() || vectEntry->first < wordEntry->first ) 
      ) {
      vector<double> valNow = vectEntry->second.valNow;
      vector<double> valDefault = vectEntry->second.valDefault;      
      if ( writeAll || valNow != valDefault ) {
        os  << vectEntry->second.name << " = ";
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
      ++vectEntry;

    // Else print word. 
    } else {
      string valNow = wordEntry->second.valNow;
      string valDefault = wordEntry->second.valDefault; 
      if ( writeAll || valNow != valDefault ) 
        os << wordEntry->second.name << " = " << valNow << "\n";
      ++wordEntry;
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
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Vect + Word Settings (all)"
       << "  -------------------------------------------------* \n";
  else if (!doListString) 
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Vect + Word Settings (cha" 
       << "nges only)  ---------------------------------------* \n" ;
  else
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Vect + Word Settings (wit" 
       << "h requested string)  ------------------------------* \n" ;
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
  map<string, Vect>::iterator vectEntry = vects.begin();
  map<string, Word>::iterator wordEntry = words.begin();

  // Loop while there is something left to do.
  while (flagEntry != flags.end() || modeEntry != modes.end() 
    || parmEntry != parms.end() || vectEntry != vects.end() 
    || wordEntry != words.end()) {

    // Check if a flag is next in lexigraphical order; if so print it.
    if ( flagEntry != flags.end() 
      && ( modeEntry == modes.end() || flagEntry->first < modeEntry->first ) 
      && ( parmEntry == parms.end() || flagEntry->first < parmEntry->first )
      && ( vectEntry == vects.end() || flagEntry->first < vectEntry->first )
      && ( wordEntry == words.end() || flagEntry->first < wordEntry->first ) 
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
      && ( vectEntry == vects.end() || modeEntry->first < vectEntry->first ) 
      && ( wordEntry == words.end() || modeEntry->first < wordEntry->first ) 
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
      && ( vectEntry == vects.end() || parmEntry->first < vectEntry->first ) 
      && ( wordEntry == words.end() || parmEntry->first < wordEntry->first ) 
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

    // Else check if vect is next, and if so print it; 
    // fixed or scientific depending on value.
    } else if ( vectEntry != vects.end()
      && ( wordEntry == words.end() || vectEntry->first < wordEntry->first ) 
      ) {
      vector<double> valsNow = vectEntry->second.valNow;
      vector<double> valsDefault = vectEntry->second.valDefault; 
      double valNow(0), valDefault(0);
      if ( doListAll || (!doListString && valsNow != valsDefault ) 
        || (doListString && vectEntry->first.find(match) != string::npos) ) {
	for (unsigned int i = 0; i < valsNow.size() || i < valsDefault.size();
	     ++i) {
	  if ( i == 0 )
	    os << " | " << setw(45) << left 
	       << vectEntry->second.name << right << " |             ";
	  else
	    os << " | " << setw(45) << " " << right << " |             ";
	  for (int j = 0; j < 4; ++j) { 
	    if (i < valsNow.size()) valNow = valsNow[i];
	    if (i < valsDefault.size()) valDefault = valsDefault[i];
	    if (j == 1) valNow = valDefault;  
	    if (j == 2) valNow = vectEntry->second.valMin;  
	    if (j == 3) valNow = vectEntry->second.valMax;  
	    if ( (j == 0 && i >= valsNow.size()) 
		 || (j == 1 && i >= valsDefault.size()) 
		 || (j == 2 && !vectEntry->second.hasMin)
		 || (j == 3 && !vectEntry->second.hasMax) )
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
      ++vectEntry;

    // Else print word. 
    } else {
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
    }
  } ;

  // End of loop over database contents.
  os << " |                                                           "
     << "                                                      | \n"
     << " *-------  End PYTHIA Flag + Mode + Parm + Vect + Word Settings  ---"
     << "-----------------------------------------------* " << endl;

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

  // Loop through the vects table, resetting all entries.
  for (map<string, Vect>::iterator vectEntry = vects.begin(); 
    vectEntry != vects.end(); ++vectEntry) {
    string name = vectEntry->first;
    resetVect(name);
  }

  // Loop through the words table, resetting all entries.
  for (map<string, Word>::iterator wordEntry = words.begin();
    wordEntry != words.end(); ++wordEntry) {
    string name = wordEntry->first;
    resetWord(name);
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

vector<double> Settings::vect(string keyIn) {
  if (isVect(keyIn)) return vects[toLower(keyIn)].valNow; 
  infoPtr->errorMsg("Error in Settings::vect: unknown key", keyIn);
  return vector<double>(1, 0.); 
}

string Settings::word(string keyIn) {
  if (isWord(keyIn)) return words[toLower(keyIn)].valNow; 
  infoPtr->errorMsg("Error in Settings::word: unknown key", keyIn);
  return " "; 
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

vector<double> Settings::vectDefault(string keyIn) {
  if (isVect(keyIn)) return vects[toLower(keyIn)].valDefault;
  infoPtr->errorMsg("Error in Settings::vectDefault: unknown key", keyIn);
  return vector<double>(1, 0.);
}

string Settings::wordDefault(string keyIn) {
  if (isWord(keyIn)) return words[toLower(keyIn)].valDefault;
  infoPtr->errorMsg("Error in Settings::wordDefault: unknown key", keyIn);
  return " ";
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

map<string, Vect> Settings::getVectMap(string match) {
  // Make the match string lower case. Start with an empty map.
  match = toLower(match);
  map<string, Vect> vectMap;  
  // Loop over the vect map (using iterator).
  for (map<string,Vect>::iterator vectEntry = vects.begin();
       vectEntry != vects.end(); ++vectEntry)
    if (vectEntry->first.find(match) != string::npos) 
      vectMap[vectEntry->first] = vectEntry->second;
  return vectMap;
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

//--------------------------------------------------------------------------
 
// Change current value, respecting limits.

void Settings::flag(string keyIn, bool nowIn) { 
  string keyLower = toLower(keyIn);
  if (isFlag(keyIn)) flags[keyLower].valNow = nowIn; 
  // Print:quiet  triggers a whole set of changes. 
  if (keyLower == "print:quiet") printQuiet( nowIn);
}

void Settings:: mode(string keyIn, int nowIn) { 
  if (isMode(keyIn)) { 
    string keyLower = toLower(keyIn);
    Mode& modeNow = modes[keyLower];
    if (modeNow.hasMin && nowIn < modeNow.valMin) 
      modeNow.valNow = modeNow.valMin; 
    else if (modeNow.hasMax && nowIn > modeNow.valMax) 
      modeNow.valNow = modeNow.valMax;
    else modeNow.valNow = nowIn; 
    // Tune:ee and Tune:pp each trigger a whole set of changes.
    if (keyLower == "tune:ee") initTuneEE( modeNow.valNow);
    if (keyLower == "tune:pp") initTunePP( modeNow.valNow); 
  } 
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

void Settings::vect(string keyIn, vector<double> nowIn) { 
  if (isVect(keyIn)) {
    Vect& vectNow = vects[toLower(keyIn)];
    vectNow.valNow.clear();
    for(vector<double>::iterator now = nowIn.begin();
	now != nowIn.end(); now++) {
      if (vectNow.hasMin && *now < vectNow.valMin) 
	vectNow.valNow.push_back(vectNow.valMin); 
      else if (vectNow.hasMax && *now > vectNow.valMax) 
	vectNow.valNow.push_back(vectNow.valMax);
      else vectNow.valNow.push_back(*now); 
    }
  } 
}  

void Settings::word(string keyIn, string nowIn) { 
    if (isWord(keyIn)) words[toLower(keyIn)].valNow = nowIn; 
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

void Settings::forceVect(string keyIn, vector<double> nowIn) { 
  if (isVect(keyIn)) vects[toLower(keyIn)].valNow = nowIn; 
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

void Settings::resetVect(string keyIn) {
  if (isVect(keyIn)) vects[toLower(keyIn)].valNow 
    = vects[toLower(keyIn)].valDefault ; 
}

void Settings::resetWord(string keyIn) {
  if (isWord(keyIn)) words[toLower(keyIn)].valNow 
    = words[toLower(keyIn)].valDefault ; 
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
  resetParm("StringZ:aLund");
  resetParm("StringZ:bLund");  
  resetParm("StringZ:rFactB");  
  resetParm("StringPT:sigma");  
  resetParm("TimeShower:alphaSvalue");  
  resetParm("TimeShower:pTmin");  
  resetParm("TimeShower:pTminChgQ"); 
 
}

//--------------------------------------------------------------------------

// Restore all pp settings to their original values.

void Settings::resetTunePP() {

  resetMode("PDF:pSet");  
  resetFlag("PDF:useLHAPDF");   
  resetParm("SigmaProcess:alphaSvalue");  
  resetFlag("SigmaTotal:zeroAXB");
  resetFlag("SigmaDiffractive:dampen");  
  resetParm("SigmaDiffractive:maxXB");
  resetParm("SigmaDiffractive:maxAX");
  resetParm("SigmaDiffractive:maxXX");  
  resetFlag("TimeShower:dampenBeamRecoil");  
  resetFlag("TimeShower:phiPolAsym");  
  resetParm("SpaceShower:alphaSvalue");  
  resetFlag("SpaceShower:samePTasMPI");  
  resetParm("SpaceShower:pT0Ref");  
  resetParm("SpaceShower:ecmRef");  
  resetParm("SpaceShower:ecmPow");  
  resetFlag("SpaceShower:rapidityOrder");  
  resetFlag("SpaceShower:phiPolAsym");  
  resetFlag("SpaceShower:phiIntAsym");  
  resetParm("MultipartonInteractions:alphaSvalue");   
  resetParm("MultipartonInteractions:pT0Ref");  
  resetParm("MultipartonInteractions:ecmRef");  
  resetParm("MultipartonInteractions:ecmPow");  
  resetMode("MultipartonInteractions:bProfile");  
  resetParm("MultipartonInteractions:expPow");  
  resetParm("MultipartonInteractions:a1");
  resetParm("BeamRemnants:primordialKTsoft");  
  resetParm("BeamRemnants:primordialKThard");  
  resetParm("BeamRemnants:halfScaleForKT");  
  resetParm("BeamRemnants:halfMassForKT");  
  resetParm("BeamRemnants:reconnectRange");
  
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
    parm("StringFlav:probStoUD",     0.30  );
    parm("StringFlav:probQQtoQ",     0.10  );
    parm("StringFlav:probSQtoQQ",    0.40  );
    parm("StringFlav:probQQ1toQQ0",  0.05  );
    parm("StringFlav:mesonUDvector", 1.00  );
    parm("StringFlav:mesonSvector",  1.50  );
    parm("StringFlav:mesonCvector",  2.50  );
    parm("StringFlav:mesonBvector",  3.00  );
    parm("StringFlav:etaSup",        1.00  );
    parm("StringFlav:etaPrimeSup",   0.40  );
    parm("StringFlav:popcornSpair",  0.50  );  
    parm("StringFlav:popcornSmeson", 0.50  );  
    parm("StringZ:aLund",            0.30  );
    parm("StringZ:bLund",            0.58  );  
    parm("StringZ:rFactB",           1.00  );  
    parm("StringPT:sigma",           0.36  );  
    parm("TimeShower:alphaSvalue",   0.137 );  
    parm("TimeShower:pTmin",         0.5   );  
    parm("TimeShower:pTminChgQ",     0.5   );  
  }

  // Marc Montull's tune to particle composition at LEP1 (August 2007).
  else if (eeTune == 2) {  
    parm("StringFlav:probStoUD",     0.22  );
    parm("StringFlav:probQQtoQ",     0.08  );
    parm("StringFlav:probSQtoQQ",    0.75  );
    parm("StringFlav:probQQ1toQQ0",  0.025 );
    parm("StringFlav:mesonUDvector", 0.5   );
    parm("StringFlav:mesonSvector",  0.6   );
    parm("StringFlav:mesonCvector",  1.5   );
    parm("StringFlav:mesonBvector",  2.5   );
    parm("StringFlav:etaSup",        0.60  );
    parm("StringFlav:etaPrimeSup",   0.15  );
    parm("StringFlav:popcornSpair",  1.0   );
    parm("StringFlav:popcornSmeson", 1.0   );
    parm("StringZ:aLund",            0.76  );
    parm("StringZ:bLund",            0.58  );   // kept fixed
    parm("StringZ:rFactB",           1.00  );   // kept fixed
    parm("StringPT:sigma",           0.36  );   // kept fixed
    parm("TimeShower:alphaSvalue",   0.137 );   // kept fixed 
    parm("TimeShower:pTmin",         0.5   );   // kept fixed 
    parm("TimeShower:pTminChgQ",     0.5   );   // kept fixed
  }

  // Full e+e- tune of flavours and FSR to LEP1 data within the 
  // Rivet + Professor framework, by Hendrik Hoeth (June 2009).
  else if (eeTune == 3) {  
    parm("StringFlav:probStoUD",     0.19  );
    parm("StringFlav:probQQtoQ",     0.09  );
    parm("StringFlav:probSQtoQQ",    1.00  );
    parm("StringFlav:probQQ1toQQ0",  0.027 );
    parm("StringFlav:mesonUDvector", 0.62  );
    parm("StringFlav:mesonSvector",  0.725 );
    parm("StringFlav:mesonCvector",  1.06  );
    parm("StringFlav:mesonBvector",  3.0   );
    parm("StringFlav:etaSup",        0.63  );
    parm("StringFlav:etaPrimeSup",   0.12  );
    parm("StringFlav:popcornSpair",  0.5   );   // kept fixed
    parm("StringFlav:popcornSmeson", 0.5   );   // kept fixed
    parm("StringZ:aLund",            0.3   );   // kept fixed
    parm("StringZ:bLund",            0.8   );  
    parm("StringZ:rFactB",           0.67  );  
    parm("StringPT:sigma",           0.304 );  
    parm("TimeShower:alphaSvalue",   0.1383);  
    parm("TimeShower:pTmin",         0.4   );   // kept fixed (near limit) 
    parm("TimeShower:pTminChgQ",     0.4   );   // kept same as pTmin
  }
  
}

//--------------------------------------------------------------------------

// Set the values related to a tune of pp/ppbar data, 
// i.e. mainly for initial-state radiation and multiparton interactions.

void Settings::initTunePP( int ppTune) {

  // Restore all pp/ppbar settings to their original values.
  // Is first step for setting up a specific tune.
  if (ppTune != 0) resetTunePP(); 

  // Decide whether to use LHAPFD where possible.
  bool preferLHAPDF = flag("Tune:preferLHAPDF");

  // Old ISR and MPI defaults from early and primitive comparisons with data.
  if (ppTune == 1) {
    mode("PDF:pSet",                            2     );  
    parm("SigmaProcess:alphaSvalue",            0.1265);  
    flag("SigmaTotal:zeroAXB",                  true  );  
    flag("SigmaDiffractive:dampen",             false );  
    flag("TimeShower:dampenBeamRecoil",         false );  
    flag("TimeShower:phiPolAsym",               false );  
    parm("SpaceShower:alphaSvalue",             0.127 );  
    flag("SpaceShower:samePTasMPI",             true  );  
    parm("SpaceShower:pT0Ref",                  2.2   );  
    parm("SpaceShower:ecmRef",                  1800.0);  
    parm("SpaceShower:ecmPow",                  0.16  );  
    flag("SpaceShower:rapidityOrder",           false );  
    flag("SpaceShower:phiPolAsym",              false );  
    flag("SpaceShower:phiIntAsym",              false );  
    parm("MultipartonInteractions:alphaSvalue", 0.127 );   
    parm("MultipartonInteractions:pT0Ref",      2.15  );  
    parm("MultipartonInteractions:ecmRef",      1800. );  
    parm("MultipartonInteractions:ecmPow",      0.16  );  
    mode("MultipartonInteractions:bProfile",    2     );  
    parm("BeamRemnants:primordialKTsoft",       0.4   );  
    parm("BeamRemnants:primordialKThard",       2.1   );  
    parm("BeamRemnants:halfScaleForKT",         7.0   );  
    parm("BeamRemnants:halfMassForKT",          2.0   );  
    parm("BeamRemnants:reconnectRange",         2.5   ); 
  }
  
  // "Tune 1" simple first tune by Peter Skands to ISR and MPI, July 2009.
  else if (ppTune == 2) {
    mode("PDF:pSet",                            2     );  
    parm("SigmaProcess:alphaSvalue",            0.1265);   
    flag("SigmaTotal:zeroAXB",                  true  );  
    flag("SigmaDiffractive:dampen",             false );  
    flag("TimeShower:dampenBeamRecoil",         false );  
    flag("TimeShower:phiPolAsym",               false );  
    parm("SpaceShower:alphaSvalue",             0.137 );  
    flag("SpaceShower:samePTasMPI",             false );  
    parm("SpaceShower:pT0Ref",                  2.0   );  
    parm("SpaceShower:ecmRef",                  1800.0);  
    parm("SpaceShower:ecmPow",                  0.0   );  
    flag("SpaceShower:rapidityOrder",           false );  
    flag("SpaceShower:phiPolAsym",              false );  
    flag("SpaceShower:phiIntAsym",              false );  
    parm("MultipartonInteractions:alphaSvalue", 0.127 );   
    parm("MultipartonInteractions:pT0Ref",      2.25  );  
    parm("MultipartonInteractions:ecmRef",      1800. );  
    parm("MultipartonInteractions:ecmPow",      0.24  );  
    mode("MultipartonInteractions:bProfile",    1     );  
    parm("BeamRemnants:primordialKTsoft",       0.5   );  
    parm("BeamRemnants:primordialKThard",       2.0   );  
    parm("BeamRemnants:halfScaleForKT",         1.0   );  
    parm("BeamRemnants:halfMassForKT",          1.0   );  
    parm("BeamRemnants:reconnectRange",         10.0  );  
  }
  
  // Tune 2C, July 2010.
  else if (ppTune == 3) {
    mode("PDF:pSet",                            8     );  
    parm("SigmaProcess:alphaSvalue",            0.135 );  
    flag("SigmaTotal:zeroAXB",                  true  );  
    flag("SigmaDiffractive:dampen",             false );  
    flag("TimeShower:dampenBeamRecoil",         true  );  
    flag("TimeShower:phiPolAsym",               true  );  
    parm("SpaceShower:alphaSvalue",             0.137 );  
    flag("SpaceShower:samePTasMPI",             false );  
    parm("SpaceShower:pT0Ref",                  2.0   );  
    parm("SpaceShower:ecmRef",                  1800.0);  
    parm("SpaceShower:ecmPow",                  0.0   );  
    flag("SpaceShower:rapidityOrder",           true  );  
    flag("SpaceShower:phiPolAsym",              true  );  
    flag("SpaceShower:phiIntAsym",              true  );  
    parm("MultipartonInteractions:alphaSvalue", 0.135 );   
    parm("MultipartonInteractions:pT0Ref",      2.32  );  
    parm("MultipartonInteractions:ecmRef",      1800. );  
    parm("MultipartonInteractions:ecmPow",      0.21  );  
    mode("MultipartonInteractions:bProfile",    3     );  
    parm("MultipartonInteractions:expPow",      1.6   );  
    parm("BeamRemnants:primordialKTsoft",       0.5   );  
    parm("BeamRemnants:primordialKThard",       2.0   );  
    parm("BeamRemnants:halfScaleForKT",         1.0   );  
    parm("BeamRemnants:halfMassForKT",          1.0   );  
    parm("BeamRemnants:reconnectRange",         3.0   );  
  }
  
  // Tune 2M, July 2010.
  else if (ppTune == 4) {
    mode("PDF:pSet",                            4     );  
    parm("SigmaProcess:alphaSvalue",            0.1265);  
    flag("SigmaTotal:zeroAXB",                  true  );  
    flag("SigmaDiffractive:dampen",             false );  
    flag("TimeShower:dampenBeamRecoil",         true  );  
    flag("TimeShower:phiPolAsym",               true  );  
    parm("SpaceShower:alphaSvalue",             0.130 );  
    flag("SpaceShower:samePTasMPI",             false );  
    parm("SpaceShower:pT0Ref",                  2.0   );  
    parm("SpaceShower:ecmRef",                  1800.0);  
    parm("SpaceShower:ecmPow",                  0.0   );  
    flag("SpaceShower:rapidityOrder",           true  );  
    flag("SpaceShower:phiPolAsym",              true  );  
    flag("SpaceShower:phiIntAsym",              true  );  
    parm("MultipartonInteractions:alphaSvalue", 0.127 );   
    parm("MultipartonInteractions:pT0Ref",      2.455 );  
    parm("MultipartonInteractions:ecmRef",      1800. );  
    parm("MultipartonInteractions:ecmPow",      0.26  );  
    mode("MultipartonInteractions:bProfile",    3     );  
    parm("MultipartonInteractions:expPow",      1.15  );  
    parm("BeamRemnants:primordialKTsoft",       0.5   );  
    parm("BeamRemnants:primordialKThard",       2.0   );  
    parm("BeamRemnants:halfScaleForKT",         1.0   );  
    parm("BeamRemnants:halfMassForKT",          1.0   );  
    parm("BeamRemnants:reconnectRange",         3.0   );  
  }
 
  // Tune 4C, October 2010.
  else if (ppTune == 5) {
    mode("PDF:pSet",                            8     );  
    parm("SigmaProcess:alphaSvalue",            0.135 );  
    flag("SigmaTotal:zeroAXB",                  true  );  
    flag("SigmaDiffractive:dampen",             true  );
    parm("SigmaDiffractive:maxXB",              65.0  );
    parm("SigmaDiffractive:maxAX",              65.0  );
    parm("SigmaDiffractive:maxXX",              65.0  );  
    flag("TimeShower:dampenBeamRecoil",         true  );  
    flag("TimeShower:phiPolAsym",               true  );  
    parm("SpaceShower:alphaSvalue",             0.137 );  
    flag("SpaceShower:samePTasMPI",             false );  
    parm("SpaceShower:pT0Ref",                  2.0   );  
    parm("SpaceShower:ecmRef",                  1800.0);  
    parm("SpaceShower:ecmPow",                  0.0   );  
    flag("SpaceShower:rapidityOrder",           true  );  
    flag("SpaceShower:phiPolAsym",              true  );  
    flag("SpaceShower:phiIntAsym",              true  );  
    parm("MultipartonInteractions:alphaSvalue", 0.135 );   
    parm("MultipartonInteractions:pT0Ref",      2.085 );  
    parm("MultipartonInteractions:ecmRef",      1800. );  
    parm("MultipartonInteractions:ecmPow",      0.19  );  
    mode("MultipartonInteractions:bProfile",    3     );  
    parm("MultipartonInteractions:expPow",      2.0   );  
    parm("BeamRemnants:primordialKTsoft",       0.5   );  
    parm("BeamRemnants:primordialKThard",       2.0   );  
    parm("BeamRemnants:halfScaleForKT",         1.0   );  
    parm("BeamRemnants:halfMassForKT",          1.0   );  
    parm("BeamRemnants:reconnectRange",         1.5   );  
  }

  // Tune 4Cx, January 2011.
  else if (ppTune == 6) {
    mode("PDF:pSet",                            8     );  
    parm("SigmaProcess:alphaSvalue",            0.135 );  
    flag("SigmaTotal:zeroAXB",                  true  );  
    flag("SigmaDiffractive:dampen",             true  );
    parm("SigmaDiffractive:maxXB",              65.0  );
    parm("SigmaDiffractive:maxAX",              65.0  );
    parm("SigmaDiffractive:maxXX",              65.0  );  
    flag("TimeShower:dampenBeamRecoil",         true  );  
    flag("TimeShower:phiPolAsym",               true  );  
    parm("SpaceShower:alphaSvalue",             0.137 );  
    flag("SpaceShower:samePTasMPI",             false );  
    parm("SpaceShower:pT0Ref",                  2.0   );  
    parm("SpaceShower:ecmRef",                  1800.0);  
    parm("SpaceShower:ecmPow",                  0.0   );  
    flag("SpaceShower:rapidityOrder",           true  );  
    flag("SpaceShower:phiPolAsym",              true  );  
    flag("SpaceShower:phiIntAsym",              true  );  
    parm("MultipartonInteractions:alphaSvalue", 0.135 );   
    parm("MultipartonInteractions:pT0Ref",      2.15  );  
    parm("MultipartonInteractions:ecmRef",      1800. );  
    parm("MultipartonInteractions:ecmPow",      0.19  );  
    mode("MultipartonInteractions:bProfile",    4     );
    parm("MultipartonInteractions:a1",          0.15  );
    parm("BeamRemnants:primordialKTsoft",       0.5   );  
    parm("BeamRemnants:primordialKThard",       2.0   );  
    parm("BeamRemnants:halfScaleForKT",         1.0   );  
    parm("BeamRemnants:halfMassForKT",          1.0   );  
    parm("BeamRemnants:reconnectRange",         1.5   );  
  }

  // Several ATLAS tunes in the A2 and AU2 series, see
  // ATLAS note ATL-PHYS-PUB-2012-003 (August 2012).
  else if (ppTune < 14) {
    parm("SigmaProcess:alphaSvalue",            0.135 );  
    flag("SigmaTotal:zeroAXB",                  true  );  
    flag("SigmaDiffractive:dampen",             true  );
    parm("SigmaDiffractive:maxXB",              65.0  );
    parm("SigmaDiffractive:maxAX",              65.0  );
    parm("SigmaDiffractive:maxXX",              65.0  );  
    flag("TimeShower:dampenBeamRecoil",         true  );  
    flag("TimeShower:phiPolAsym",               true  );  
    parm("SpaceShower:alphaSvalue",             0.137 );  
    flag("SpaceShower:samePTasMPI",             false );  
    parm("SpaceShower:pT0Ref",                  2.0   );  
    parm("SpaceShower:ecmRef",                  1800.0);  
    parm("SpaceShower:ecmPow",                  0.0   );  
    flag("SpaceShower:rapidityOrder",           false );  
    flag("SpaceShower:phiPolAsym",              true  );  
    flag("SpaceShower:phiIntAsym",              true  );  
    parm("MultipartonInteractions:alphaSvalue", 0.135 );   
    parm("MultipartonInteractions:ecmRef",      1800. );  
    mode("MultipartonInteractions:bProfile",    4     );
    parm("BeamRemnants:primordialKTsoft",       0.5   );  
    parm("BeamRemnants:primordialKThard",       2.0   );  
    parm("BeamRemnants:halfScaleForKT",         1.0   );  
    parm("BeamRemnants:halfMassForKT",          1.0   ); 

    // ATLAS MB tune A2-CTEQ6L1.
    if (ppTune == 7) { 
      if (preferLHAPDF) {
        flag("PDF:useLHAPDF",                   true  );   
        word("PDF:LHAPDFset",          "cteq6ll.LHpdf");
      } else mode("PDF:pSet",                   8     );  
      parm("MultipartonInteractions:pT0Ref",    2.18  );  
      parm("MultipartonInteractions:ecmPow",    0.22  );  
      parm("MultipartonInteractions:a1",        0.06  );
      parm("BeamRemnants:reconnectRange",       1.55  );  
    }

    // ATLAS MB tune A2-MSTW2008LO.
    else if (ppTune == 8) { 
      if (preferLHAPDF) {
        flag("PDF:useLHAPDF",                   true  );   
        word("PDF:LHAPDFset",  "MSTW2008lo68cl.LHgrid");
      } else mode("PDF:pSet",                   5     );  
      parm("MultipartonInteractions:pT0Ref",    1.90  );  
      parm("MultipartonInteractions:ecmPow",    0.30  );  
      parm("MultipartonInteractions:a1",        0.03  );
      parm("BeamRemnants:reconnectRange",       2.28  );  
    }

    // ATLAS UE tune AU2-CTEQ6L1.
    if (ppTune == 9) { 
      if (preferLHAPDF) {
        flag("PDF:useLHAPDF",                   true  );   
        word("PDF:LHAPDFset",          "cteq6ll.LHpdf");
      } else mode("PDF:pSet",                   8     );  
      parm("MultipartonInteractions:pT0Ref",    2.13  );  
      parm("MultipartonInteractions:ecmPow",    0.21  );  
      parm("MultipartonInteractions:a1",        0.00  );
      parm("BeamRemnants:reconnectRange",       2.21  );  
    }

    // ATLAS UE tune AU2-MSTW2008LO.
    else if (ppTune == 10) { 
      if (preferLHAPDF) {
        flag("PDF:useLHAPDF",                   true  );   
        word("PDF:LHAPDFset",  "MSTW2008lo68cl.LHgrid");
      } else mode("PDF:pSet",                   5     );  
      parm("MultipartonInteractions:pT0Ref",    1.87  );  
      parm("MultipartonInteractions:ecmPow",    0.28  );  
      parm("MultipartonInteractions:a1",        0.01  );
      parm("BeamRemnants:reconnectRange",       5.32  );  
    }

    // ATLAS UE tune AU2-CT10.
    else if (ppTune == 11) { 
      flag("PDF:useLHAPDF",                     true  );   
      word("PDF:LHAPDFset",              "CT10.LHgrid");
      parm("MultipartonInteractions:pT0Ref",    1.70  );  
      parm("MultipartonInteractions:ecmPow",    0.16  );  
      parm("MultipartonInteractions:a1",        0.10  );
      parm("BeamRemnants:reconnectRange",       4.67  );  
    }

    // ATLAS UE tune AU2-MRST2007LO*.
    else if (ppTune == 12) { 
      if (preferLHAPDF) {
        flag("PDF:useLHAPDF",                   true  );   
        word("PDF:LHAPDFset",   "MRST2007lomod.LHgrid");
      } else mode("PDF:pSet",                   3     );  
      parm("MultipartonInteractions:pT0Ref",    2.39  );  
      parm("MultipartonInteractions:ecmPow",    0.24  );  
      parm("MultipartonInteractions:a1",        0.01  );
      parm("BeamRemnants:reconnectRange",       1.76  );  
    }

    // ATLAS UE tune AU2-MRST2007LO**.
    else if (ppTune == 13) { 
      if (preferLHAPDF) {
        flag("PDF:useLHAPDF",                   true  );   
        word("PDF:LHAPDFset",        "MRSTMCal.LHgrid");
      } else mode("PDF:pSet",                   4     );  
      parm("MultipartonInteractions:pT0Ref",    2.57  );  
      parm("MultipartonInteractions:ecmPow",    0.23  );  
      parm("MultipartonInteractions:a1",        0.01  );
      parm("BeamRemnants:reconnectRange",       1.47  );  
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

// Extract XML vector value following XML attribute.

vector<double> Settings::vectorAttributeValue(string line, string attribute) {
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
