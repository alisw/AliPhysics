//------------------------------------------------------------------------------
// Implementation of abstract AliComparisonObject class. It keeps information from 
// comparison of reconstructed and MC particle tracks. All comparison objects should 
// derive from this class.
//
// Author: J.Otwinowski 14/04/2008 
//------------------------------------------------------------------------------

#include <iostream>

#include "AliMCInfo.h" 
#include "AliESDRecInfo.h" 
#include "AliComparisonObject.h" 

using namespace std;

ClassImp(AliComparisonObject)

//_____________________________________________________________________________
AliComparisonObject::AliComparisonObject():
  TNamed("AliComparisonObject","AliComparisonObject") {
  // constructor
}

//_____________________________________________________________________________
AliComparisonObject::AliComparisonObject(const char* name):
  TNamed(name,name) {
  // constructor
}

//_____________________________________________________________________________
AliComparisonObject::~AliComparisonObject(){
  // destructor 
}

