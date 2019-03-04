#include "AliAnalysisTaskMKTest.h"

class AliAnalysisTaskMKTest;

ClassImp(AliAnalysisTaskMKTest)

//_____________________________________________________________________________

AliAnalysisTaskMKTest::AliAnalysisTaskMKTest() 
    : AliAnalysisTaskSE()     
{
    //constructor
}

//_____________________________________________________________________________

AliAnalysisTaskMKTest::AliAnalysisTaskMKTest(const char* name)
    : AliAnalysisTaskSE(name)
{
    // constructor
//      DefineInput(0, TChain::Class());    // define input
//      DefineOutput(1, TList::Class());    // define ouptut 
}

//_____________________________________________________________________________

AliAnalysisTaskMKTest::~AliAnalysisTaskMKTest()
{
    // destructor
   
}
