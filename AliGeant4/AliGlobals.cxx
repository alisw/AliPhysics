// $Id$
// Category: global
//
// See the class description in the header file.

#include "AliGlobals.h"
#include "TG4Globals.h"

#include <stdlib.h>

// static data members
const G4double AliGlobals::fgkDefaultCut = 2.0*mm;

AliGlobals::AliGlobals() {
//
}
  
AliGlobals::~AliGlobals() {
//
}
  
// static methods

void AliGlobals::Exception(const char* s)
{
// Prints error message end exits the program.
// ---

  if (s)
  {  G4cerr << G4endl << "    " << s << G4endl; }
  G4cerr << "*** AliceException: Aborting execution ***" << G4endl;   
  exit(1);
}

void AliGlobals::Warning(const char* s)
{
// Prints warning message.
// ---

  G4cerr << "+++ Alice Warning: +++" << G4endl;   
  if (s)
  {  G4cerr  << "    " << s << G4endl; }
  G4cerr << "++++++++++++++++++++++" << G4endl;   
}

#ifdef G4USE_STL
void AliGlobals::Exception(G4std::string s) {
//
  AliGlobals::Exception(s.c_str());
}

void AliGlobals::Exception(G4String s) {
//
   AliGlobals::Exception(s.c_str());
}

void AliGlobals::Warning(G4std::string s) {
//
  AliGlobals::Warning(s.c_str());
}

void AliGlobals::Warning(G4String s) {
//
  AliGlobals::Warning(s.c_str());
}
#endif

void AliGlobals::AppendNumberToString(G4String& s, G4int a)
{
// Appends number to string.
// ---
  TG4Globals::AppendNumberToString(s, a);
}

G4int AliGlobals::StringToInt(G4String s)
{
// Converts one char string to integer number.
// ---

  // make better
  if (s=="0") return 0;
  if (s=="1") return 1;
  if (s=="2") return 2;
  if (s=="3") return 3;
  if (s=="4") return 4;
  if (s=="5") return 5;
  if (s=="6") return 6;
  if (s=="7") return 7;
  if (s=="8") return 8;
  if (s=="9") return 9;
  return -1;
}
 
