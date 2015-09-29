/**********************************************
 *
 * Class designed to store all multiplicity
 * variables in a TClonesArray
 * 
 * Instancing an empty AliMultInput class
 * will allow you to store standard variables
 *
 **********************************************/

#include "TList.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"

ClassImp(AliMultInput);

AliMultInput::AliMultInput() :
  TNamed(), fNVars(0), fVariableList(0x0)
{
  // Constructor
    fVariableList = new TList();
}
AliMultInput::AliMultInput(const char * name, const char * title):
TNamed(name,title), fNVars(0), fVariableList(0x0)
{
  // Constructor
    fVariableList = new TList();
}
AliMultInput::~AliMultInput(){
  // destructor
  
}