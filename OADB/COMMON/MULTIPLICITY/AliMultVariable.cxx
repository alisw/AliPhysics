/**********************************************
 *
 * This is the basic input variable definition 
 * class
 * 
 * Requires a definition of a callsign and 
 * a definition of a substitution value
 *
 **********************************************/

#include "AliMultVariable.h"

ClassImp(AliMultVariable);

AliMultVariable::AliMultVariable() :
  TNamed(), fIsInteger(kFALSE), fValue(0), fValueInteger(0), fMean(0)
{
  // Constructor
  
}
AliMultVariable::AliMultVariable(const char * name, const char * title):
TNamed(name,title), fIsInteger(kFALSE), fValue(0), fValueInteger(0), fMean(0)
{
  // Constructor
  
}
AliMultVariable::~AliMultVariable(){
  // destructor
  
}