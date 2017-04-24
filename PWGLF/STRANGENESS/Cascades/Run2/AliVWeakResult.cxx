//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold weak decay results
// This is a base class for AliVWeakResult and AliCascadeResult
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include "TList.h"
#include "TH3F.h"
#include "TProfile.h"
#include "AliVWeakResult.h"
#include "AliLog.h"
#include <iostream>
#include <TROOT.h>
using namespace std;

ClassImp(AliVWeakResult);

//________________________________________________________________
AliVWeakResult::AliVWeakResult() :
TNamed()
{
    // Dummy Constructor - not to be used!
}

//________________________________________________________________
AliVWeakResult::AliVWeakResult(const char * name, const char * title) :
TNamed(name,title)
{
    // TNamed-inspired constructor
}

//________________________________________________________________
AliVWeakResult::~AliVWeakResult(){
    // A boring class. Nothing to delete.
}

