//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold Machines
// virtual class: other specific classes will derive from this
// but share the interface!
// This is a base class for AliNeuralNetwork and AliBDT
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include "AliMachineLearning.h"
#include <iostream>
#include <TROOT.h>
using namespace std;

ClassImp(AliMachineLearning);

//________________________________________________________________
AliMachineLearning::AliMachineLearning() : TNamed()
{
    // Dummy Constructor - not to be used!
}

//________________________________________________________________
AliMachineLearning::AliMachineLearning(const char * name, const char * title) : TNamed(name,title)
{
    // TNamed-inspired constructor
}

//________________________________________________________________
AliMachineLearning::~AliMachineLearning()
{
    // A boring class. Nothing to delete.
}

