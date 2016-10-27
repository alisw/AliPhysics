/**********************************************
 *
 * Class designed to provide a basic published
 * interface for AliRoot classes to be able 
 * to reference centrality information
 * 
 * For completeness, the actual methods that
 * will typically be invoked at runtime are the 
 * ones in the class: 
 * 
 * OADB/COMMON/MULTIPLICITY/AliMultSelection.cxx
 * 
 * First implementation includes:
 * --- GetMultiplicityPercentile
 * --- GetEvSelCode 
 *
 * Bugs? Problems? Suggestions? Please contact:
 * --- david.dobrigkeit.chinellato@cern.ch 
 *
 **********************************************/

#include <iostream>
#include <TROOT.h>
#include "AliMultSelectionBase.h"
using namespace std;

ClassImp(AliMultSelectionBase);
//________________________________________________________________
AliMultSelectionBase::AliMultSelectionBase() :
  TNamed()
{
  // Constructor
}
//________________________________________________________________
AliMultSelectionBase::AliMultSelectionBase(const char * name, const char * title):
TNamed(name,title)
{
  // Constructor
}
//________________________________________________________________
AliMultSelectionBase::~AliMultSelectionBase(){
    // destructor: clean stuff up
    //Nothing to destroy
}