// $Id$
// Category: global
//
// See the class description in the header file.

#include "TG4Globals.h"

#include <stdlib.h>

TG4Globals::TG4Globals() {
//
}
  
TG4Globals::~TG4Globals() {
//
}
  
// static methods

void TG4Globals::Exception(const char* string)
{
// Prints error message end exits the program.
// ---

  if (string)
  {  G4cerr << G4endl << "    " << string << G4endl; }
  G4cerr << "*** TG4Exception: Aborting execution ***" << G4endl;   
  exit(1);
}

void TG4Globals::Warning(const char* string)
{
// Prints warning message.
// ---

  G4cerr << "++++  TG4Warning:  ++++" << G4endl;   
  if (string)
  {  G4cerr  << "    " << string << G4endl; }
  G4cerr << "+++++++++++++++++++++++" << G4endl;   
}

void TG4Globals::AppendNumberToString(G4String& s, G4int a)
{
// Appends number to string.
// ---

  const char* kpNumber="0123456789";
  G4String p=""; G4String q="";
  do 
  {
    G4int b=a/10;
    G4int c=a%10;
    p=kpNumber[c];
    q=p.append(q);
    a=b;        
  } while (a>0);
  s.append(q);
}

