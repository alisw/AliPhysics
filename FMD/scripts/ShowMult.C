//
// $Id$
//
// Script to dump multiplicity information to std::cout. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <AliFMDMultStrip.h>
#include <AliFMDMultRegion.h>
#include <AliFMDInput.h>

class ShowMult : public AliFMDInputRecPoints
{
public:
  ShowMult() {}
  Bool_t ProcessStrip(AliFMDMultStrip* mult) 
  {
    mult->Print();
    return kTRUE;
  }
  Bool_t ProcessRegion(AliFMDMultRegion* mult) 
  {
    mult->Print();
    return kTRUE;
  }
};

//____________________________________________________________________
//
// EOF
//
