// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4Globals
// ----------------
// See the class description in the header file.

#include "TG4Globals.h"

#include <stdlib.h>

//_____________________________________________________________________________
TG4Globals::TG4Globals() {
//
}
  
//_____________________________________________________________________________
TG4Globals::~TG4Globals() {
//
}
  
// static methods

//_____________________________________________________________________________
void TG4Globals::Exception(const char* string)
{
// Prints error message end exits the program.
// ---

  if (string)
  {  G4cerr << G4endl << "    " << string << G4endl; }
  G4cerr << "*** TG4Exception: Aborting execution ***" << G4endl;   
  exit(1);
}

//_____________________________________________________________________________
void TG4Globals::Warning(const char* string)
{
// Prints warning message.
// ---

  G4cerr << "++++  TG4Warning:  ++++" << G4endl;   
  if (string)
  {  G4cerr  << "    " << string << G4endl; }
  G4cerr << "+++++++++++++++++++++++" << G4endl;   
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
G4bool TG4Globals::Compare(G4bool activation, TG4G3ControlValue controlValue)
{
// Compares the boolean value of the process activation
// with the process control value.
// Returns true if the values correspond, false otherwise.
// ---

  if (controlValue == kUnset) {
    TG4Globals::Warning(
      "TG4SpecialControls::Compare: control value = kUnset.");
    return false;
  }    

  if (controlValue == kActivate || controlValue == kActivate2)
    return activation;
  else
    return !activation;  
}  

//_____________________________________________________________________________
void TG4Globals::PrintStars(G4bool emptyLineFirst)
{
// Print stars.
// ---
  

  if (emptyLineFirst)  G4cout << G4endl;
  
  G4cout << "**********************************************" << G4endl;
     
  if (!emptyLineFirst) G4cout << G4endl;
}

