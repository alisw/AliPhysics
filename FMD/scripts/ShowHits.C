//
// $Id$
//
// Script to dump hit information to std::cout. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <TH2D.h>
#include <AliFMDHit.h>
#include <AliFMDInput.h>

class ShowHits : public AliFMDInputHits
{
public:
  ShowHits() {}
  Bool_t ProcessHit(AliFMDHit* hit, TParticle* part) 
  {
    if (!hit) {
      std::cout << "No hit" << std::endl;
      return kFALSE;
    }
    hit->Print();
    if (part) part->Print();
    return kTRUE;
  }
};

//____________________________________________________________________
//
// EOF
//
