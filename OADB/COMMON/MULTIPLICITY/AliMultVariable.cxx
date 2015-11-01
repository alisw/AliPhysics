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

//________________________________________________________________
AliMultVariable::AliMultVariable() :
  TNamed(), fIsInteger(kFALSE), fValue(0), fValueInteger(0), fMean(0)
{
  // Constructor
  
}
//________________________________________________________________
AliMultVariable::AliMultVariable(const char * name, const char * title):
TNamed(name,title), fIsInteger(kFALSE), fValue(0), fValueInteger(0), fMean(0)
{
  // Constructor
  
}
//________________________________________________________________
AliMultVariable::~AliMultVariable(){
  // destructor
  
}
//________________________________________________________________
void AliMultVariable::Print(Option_t* option) const
{
    printf("%s: %s/%s (%s)",
           ClassName(), GetName(), GetTitle(),
           (IsInteger() ? "integer" : "float"));
    TString o(option);
    if (o.Contains("V")){
        if (IsInteger()) printf("=%d", fValueInteger);
        else             printf("=%f", fValue);
    }
    if (o.Contains("M")) printf(" (mean=%f)", fMean);
    printf("\n");
}