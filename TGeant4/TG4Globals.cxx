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
  {  cerr << endl << "    " << string << endl; }
  cerr << "*** TG4Exception: Aborting execution ***" << endl;   
  exit(1);
}

void TG4Globals::Warning(const char* string)
{
// Prints warning message.
// ---

  cerr << "++++  TG4Warning:  ++++" << endl;   
  if (string)
  {  cerr  << "    " << string << endl; }
  cerr << "+++++++++++++++++++++++" << endl;   
}
