#include "AliAnalysis.h"
//________________________________
///////////////////////////////////////////////////////////
//
// class AliAnalysis
//
// Base class for analysis.
// Each inheriting calss must define 3 methods:
//   - Init() : that is called before event processing
//   - ProcessEvent(AliESD*,AliStack*)
//   - 
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////


ClassImp(AliAnalysis)

AliAnalysis::AliAnalysis()
{
 //ctor
}
/*********************************************************/

AliAnalysis::AliAnalysis(const char* name,const char* title):
 TTask(name,title)
{
 //ctor
}
/*********************************************************/

AliAnalysis::~AliAnalysis()
{
 //dtor
}
/*********************************************************/

