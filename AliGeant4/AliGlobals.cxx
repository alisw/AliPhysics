// $Id$
// Category: global
//
// See the class description in the header file.

#include "AliGlobals.h"

#include <stdlib.h>

// static data members
const G4double AliGlobals::fgDefaultCut = 2.0*mm;

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
  {  cerr << endl << "    " << s << endl; }
  cerr << "*** AliceException: Aborting execution ***" << endl;   
  exit(1);
}

void AliGlobals::Warning(const char* s)
{
// Prints warning message.
// ---

  cerr << "+++ Alice Warning: +++" << endl;   
  if (s)
  {  cerr  << "    " << s << endl; }
  cerr << "++++++++++++++++++++++" << endl;   
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
};

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
 
