#include "AliAnalysis.h"
//________________________________
///////////////////////////////////////////////////////////
//
// class AliAnalysis
//
// Base class for analysis
//
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////


ClassImp(AliAnalysis)

Int_t AliAnalysis::fgkDebug = 0;

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

